##
## Show progress and estimated time of arrival to the stderr.
## Requires to know the number of units in advance.
## Prints a message at given steps of fraction of units completed.
##
import strformat
import times

type
  Progress* = object
    n_units: int
    n_completed: int
    f_completed: float
    last_printed_f: float
    step: float
    start_time: Time
    outfile: File
    info: string

proc new_progress*(n_units: int, step: float = 0.01,
                   info: string = ""): Progress =
  result.n_units = n_units
  result.step = step
  result.start_time = get_time()
  result.outfile = stderr
  result.info = info

proc `$`*(self: Progress): string =
  let
    elapsed = get_time()-self.start_time
    total = ((elapsed.in_milliseconds.float/self.f_completed)/1000).int
    toa = total - elapsed.in_seconds
    infostr = if len(self.info) > 0: self.info & "; " else: ""
  result =
    &"# {infostr}completed {self.n_completed} of {self.n_units} "&
    &"({(self.f_completed*100):.1f}% completed; " &
    &"elapsed: {elapsed.in_seconds}s; TOA: {toa}s; total: {total}s)"

proc inc*(self: var Progress, quiet = false) =
  ## increments the counter of completed units and prints the information
  ## to the stderr when a new step is completed
  self.n_completed += 1
  self.f_completed = self.n_completed/self.n_units
  if not quiet and self.f_completed == 1.0 or
      self.f_completed >= self.last_printed_f + self.step:
    self.outfile.write_line($self)
    self.last_printed_f = self.f_completed
