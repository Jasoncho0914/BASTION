# Calendar preprocessing for data already indexed by POSIXct dates.
# Parsing text dates is the caller's responsibility. No fitting, aggregation,
# interpolation, or guessing of the data collection interval occurs here.

.calendar_dates <- function(dates,tz=NULL,label="dates") {
  if(!inherits(dates,"POSIXct") || !is.null(dim(dates)) || !length(dates))
    stop(label," must be a nonempty POSIXct vector; convert text dates before calling")
  if(any(!is.finite(as.numeric(dates))))
    stop(label," contains missing/invalid timestamp(s) at row(s): ",paste(which(!is.finite(as.numeric(dates))),collapse=", "))
  if(is.null(tz)) {
    z<-attr(dates,"tzone")
    if(length(z) && !is.na(z[1]) && nzchar(z[1]))tz<-z[1]
  }
  if(is.null(tz) || length(tz)!=1L || is.na(tz) || !(tz%in%c("UTC","GMT",OlsonNames())))
    stop("Set a valid time zone on ",label," or supply tz; no system-local zone is assumed")
  list(time=lubridate::with_tz(dates,tz),tz=tz)
}

.calendar_seconds <- function(spec) {
  spec$every*switch(spec$unit,millisecond=.001,second=1,minute=60,hour=3600)
}

.calendar_step <- function(spec,label,seasonal=FALSE) {
  # Optional positional shorthand: list("daily", every=1, ...).
  if(is.list(spec) && length(spec) && (is.null(names(spec)) || is.na(names(spec)[1]) || names(spec)[1]=="")) {
    if(!is.character(spec[[1]]) || length(spec[[1]])!=1L || !is.null(spec$unit))
      stop(label," positional first entry must be a unit string, without another unit field")
    first<-spec[[1]];spec<-spec[-1];spec$unit<-first
  }
  if(!is.list(spec) || is.null(spec$unit))stop(label," must be a list containing unit and optional every/anchor")
  extra<-setdiff(names(spec),c("unit","every","anchor",if(seasonal)"distance_scale"))
  if(length(extra))stop("Unknown ",label," field(s): ",paste(extra,collapse=", "))
  unit<-spec$unit
  aliases<-c(daily="day",weekly="week",monthly="month",quarterly="quarter",yearly="year",hourly="hour",minutely="minute",secondly="second")
  if(is.character(unit) && length(unit)==1L && !is.na(unit) && unit%in%names(aliases))unit<-unname(aliases[unit])
  allowed<-c("millisecond","second","minute","hour","day","week","month","quarter","year")
  if(!is.character(unit) || length(unit)!=1L || is.na(unit) || !(unit%in%allowed))
    stop(label," unit must be one of: ",paste(allowed,collapse=", "))
  every<-if(is.null(spec$every))1 else spec$every
  if(!is.numeric(every) || length(every)!=1L || !is.finite(every) || every<=0)
    stop(label," every must be a positive finite number")
  calendar<-unit%in%c("day","week","month","quarter","year")
  if(calendar && every!=round(every))stop("Calendar every must be a positive integer")
  result<-list(unit=unit,every=every,anchor=spec$anchor,calendar=calendar)
  if(seasonal) {
    scaled<-if(is.null(spec$distance_scale))TRUE else spec$distance_scale
    if(!is.logical(scaled) || length(scaled)!=1L || is.na(scaled))stop(label," distance_scale must be TRUE or FALSE")
    result$distance_scale<-scaled
  }
  result
}

# Calendar arithmetic is from the ORIGINAL anchor, never repeated month-end
# rollovers. An impossible requested date is an error, not a clamped date.
.calendar_clock <- function(spec,anchor,tz) {
  unit<-spec$unit;every<-spec$every
  at<-if(!spec$calendar){
    seconds<-every*switch(unit,millisecond=.001,second=1,minute=60,hour=3600)
    function(n)as.POSIXct(as.numeric(anchor)+n*seconds,origin="1970-01-01",tz=tz)
  } else {
    amount<-every*switch(unit,day=1,week=7,month=1,quarter=3,year=12)
    period_unit<-if(unit%in%c("day","week"))"days" else "months"
    function(n){
      ans<-suppressWarnings(anchor+(if(period_unit=="days")lubridate::days(n*amount) else lubridate::period(month=n*amount)))
      if(anyNA(ans))stop("Anchor ",format(anchor,tz=tz)," produces an impossible calendar date; choose another anchor (no automatic month-end rollover)")
      ans
    }
  }
  locate<-function(times){
    if(!spec$calendar)q<-floor((as.numeric(times)-as.numeric(anchor))/seconds)
    else if(unit%in%c("day","week"))q<-floor(as.numeric(as.Date(times,tz=tz)-as.Date(anchor,tz=tz))/amount)
    else {
      month_index<-function(x)12*lubridate::year(x)+lubridate::month(x)
      q<-floor((month_index(times)-month_index(anchor))/amount)
    }
    repeat {
      left<-as.numeric(at(q));right<-as.numeric(at(q+1))
      down<-as.numeric(times)<left;up<-as.numeric(times)>=right
      if(!any(down|up))break
      q[down]<-q[down]-1;q[up]<-q[up]+1
    }
    q
  }
  list(at=at,locate=locate,anchor=anchor,spec=spec)
}

.calendar_anchor <- function(spec,parsed,observation=FALSE) {
  if(!is.null(spec$anchor)) {
    if(length(spec$anchor)!=1L)stop("anchor must have length one")
    return(.calendar_dates(spec$anchor,tz=parsed$tz,label="anchor")$time)
  }
  if(observation)return(min(parsed$time))
  # Reference origin for multi-unit cycles; Monday for weekly cycles.
  value<-if(spec$unit=="week")"1970-01-05" else "1970-01-01"
  as.POSIXct(value,tz=parsed$tz)
}

# One data vector and one POSIXct vector; every seasonal component gets its
# own calendar layout. Both data_unit and seasonalities use unit/every/anchor.
# Calendar defaults:
# every=1; seasonal reference Jan 1 1970 at local midnight (weeks: Monday
# Jan 5); data collection anchor=earliest supplied date. All returned explicitly.
# tolerance_seconds=NULL uses 2*machine epsilon*max(1,abs(dates),abs(data anchor)).
# Anchors, when supplied, must also be POSIXct. No character parsing interface.
#' @noRd
prepare_calendar <- function(data,dates,data_unit,seasonalities,tz=NULL,
                             tolerance_seconds=NULL) {
  if(!requireNamespace("lubridate",quietly=TRUE))stop("Install lubridate to use calendar preprocessing")
  parsed<-.calendar_dates(dates,tz)
  if(!is.numeric(data) || !is.null(dim(data)) || length(data)!=length(dates) || any(is.infinite(data)))
    stop("data must be a numeric vector with one finite value or NA per date")
  if(!is.list(seasonalities) || !length(seasonalities) || !all(vapply(seasonalities,is.list,logical(1))))
    stop("seasonalities must be a nonempty list of seasonal specification lists")
  obs<-.calendar_step(data_unit,"data_unit")
  obs_anchor<-.calendar_anchor(obs,parsed,observation=TRUE)
  obs_clock<-.calendar_clock(obs,obs_anchor,parsed$tz)
  if(!is.null(tolerance_seconds) && (!is.numeric(tolerance_seconds) || length(tolerance_seconds)!=1L ||
      !is.finite(tolerance_seconds) || tolerance_seconds<0))stop("tolerance_seconds must be NULL or a nonnegative finite scalar")
  tol<-if(is.null(tolerance_seconds))2*.Machine$double.eps*max(1,abs(as.numeric(parsed$time)),abs(as.numeric(obs_anchor))) else tolerance_seconds
  if(!obs$calendar && .calendar_seconds(obs)<=2*tol)
    stop("data_unit interval is too small for POSIXct timestamp precision or tolerance_seconds")
  labels<-names(seasonalities)
  if(is.null(labels))labels<-rep("",length(seasonalities))
  layouts<-vector("list",length(seasonalities))
  for(p in seq_along(seasonalities)){
    spec<-.calendar_step(seasonalities[[p]],paste0("seasonalities[[",p,"]]"),seasonal=TRUE)
    if(!spec$calendar && .calendar_seconds(spec)<=2*tol)
      stop("Seasonal interval is too small for POSIXct timestamp precision or tolerance_seconds")
    anchor<-.calendar_anchor(spec,parsed)
    clock<-.calendar_clock(spec,anchor,parsed$tz)
    cycles<-clock$locate(parsed$time);indices<-seq(min(cycles),max(cycles))
    boundaries<-clock$at(c(indices,max(indices)+1))
    if(any(diff(as.numeric(boundaries))<=0))stop("Cycle boundaries must be strictly increasing")
    # Construct observation slots covering every complete cycle, including the
    # unobserved portions of first/last cycles and wholly unobserved cycles.
    ends<-obs_clock$locate(boundaries[c(1,length(boundaries))])
    slots<-obs_clock$at(seq(min(ends)-1,max(ends)+1))
    slots<-slots[as.numeric(slots)>=as.numeric(boundaries[1]) & as.numeric(slots)<as.numeric(tail(boundaries,1))]
    if(!length(slots))stop("No observation slots within requested cycles")
    cycle_id<-findInterval(as.numeric(slots),as.numeric(boundaries))
    lengths<-tabulate(cycle_id,nbins=length(indices))
    if(any(lengths<2L))stop("Requested seasonality is too short for the declared data collection interval: fewer than two slots in a cycle")
    if(length(slots)>1L && tol>=min(diff(as.numeric(slots)))/2)
      stop("Timestamp tolerance is too large relative to observation spacing")
    # Match to nearest slot, checking displacement explicitly.
    left<-findInterval(as.numeric(parsed$time),as.numeric(slots))
    a<-pmax(1L,pmin(length(slots),left));b<-pmax(1L,pmin(length(slots),left+1L))
    use_b<-abs(as.numeric(parsed$time)-as.numeric(slots[b]))<abs(as.numeric(parsed$time)-as.numeric(slots[a]))
    matched<-a;matched[use_b]<-b[use_b]
    displacement<-abs(as.numeric(parsed$time)-as.numeric(slots[matched]))
    if(any(displacement>tol))stop("Timestamp(s) off the declared data collection interval at original row(s): ",paste(which(displacement>tol),collapse=", "))
    if(anyDuplicated(matched))stop("Multiple observations occupy the same declared slot; no automatic aggregation is performed")
    within<-sequence(lengths)
    phase<-(within-.5)/lengths[cycle_id]
    grid<-sort(unique(phase));grid_index<-match(phase,grid)
    if(length(grid)<3L)stop("Current second-difference seasonal sampler requires at least three distinct latent phase positions")
    observed_positions<-unique(grid_index[matched])
    if(length(observed_positions)<2L)stop("Observed timestamps cover only one seasonal phase; the requested seasonal shape is not identifiable from this schedule")
    original_row<-rep(NA_integer_,length(slots));original_row[matched]<-seq_along(parsed$time)
    label<-labels[p];if(is.na(label) || !nzchar(label))label<-paste0(spec$unit,"_",spec$every)
    layout<-list(name=label,specification=spec,distance_scale=spec$distance_scale,resolved_anchor=anchor,tz=parsed$tz,
      boundaries=boundaries,absolute_cycle_index=indices,grid=grid,cycle_id=cycle_id,
      grid_index=grid_index,cycle_lengths=lengths,slot_time=slots,
      original_to_slot=matched,slot_to_original=original_row,
      slots=data.frame(slot=seq_along(slots),time=slots,cycle_id=cycle_id,
        position=within,cycle_length=lengths[cycle_id],phase=phase,
        grid_index=grid_index,original_row=original_row),
      n_original=length(parsed$time),n_slots=length(slots),n_missing_slots=sum(is.na(original_row)),
      grid_size=length(grid),timestamp_displacement_seconds=displacement)
    layout$data<-calendar_expand(data,layout)
    layouts[[p]]<-layout
  }
  names(layouts)<-vapply(layouts,`[[`,character(1),"name")
  if(anyDuplicated(names(layouts)))stop("Seasonality names must be unique; name components with different anchors explicitly")
  structure(list(data=as.numeric(data),dates=parsed$time,tz=parsed$tz,
    data_unit=obs,data_anchor=obs_anchor,
    tolerance_seconds=tol,seasonalities=layouts),class="bastion_calendar")
}

# Insert explicit missing values without modifying the supplied observations.
#' @noRd
calendar_expand <- function(values,layout) {
  if(!is.numeric(values) || length(values)!=layout$n_original || any(is.infinite(values)))
    stop("values must be numeric with one finite value or NA per original timestamp")
  ans<-rep(NA_real_,layout$n_slots);ans[layout$original_to_slot]<-values;ans
}

# Return evaluated seasonal estimates to the user's original observation order.
#' @noRd
calendar_restore <- function(values,layout) {
  if(is.null(dim(values))) {
    if(length(values)!=layout$n_slots)stop("Expected one value per complete-layout slot")
    return(values[layout$original_to_slot])
  }
  if(nrow(values)!=layout$n_slots)stop("Matrix rows must correspond to complete-layout slots")
  values[layout$original_to_slot,,drop=FALSE]
}

#' @noRd
print.bastion_calendar <- function(x,...) {
  cat("Calendar preprocessing\nTime zone:",x$tz,"\nData collection interval:",
    x$data_unit$every,x$data_unit$unit,"\nData collection anchor:",
    base::format(x$data_anchor,tz=x$tz),"\nTimestamp tolerance (seconds):",x$tolerance_seconds,"\n")
  for(d in x$seasonalities)cat("\n",d$name,": ",length(d$cycle_lengths)," cycles; ",d$n_slots,
    " slots; ",d$n_missing_slots," missing slots; grid size ",d$grid_size,
    "\n  cycle lengths: ",paste(sort(unique(d$cycle_lengths)),collapse=", "),
    "; anchor: ",base::format(d$resolved_anchor,tz=d$tz),"\n",sep="")
  invisible(x)
}
