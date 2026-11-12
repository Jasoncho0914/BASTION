# Shared progress display for fixed-period and calendar fits.
.bastion_progress_bar <- function(total) {
  progress::progress_bar$new(
    format="(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total=total,complete="=",incomplete="-",current=">",clear=FALSE,width=100)
}

.bastion_progress_tick <- function(bar,iteration,total) {
  if(iteration<10L || iteration%%100L==0L || iteration==total)
    bar$update(iteration/total)
  invisible(NULL)
}

.bastion_progress_chain <- function(chain,verbose) {
  if(verbose)print(paste("Chain",chain))
  invisible(NULL)
}
