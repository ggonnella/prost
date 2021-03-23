import macros
import nimpy

macro exportconst(cnst: untyped): untyped =
  let procname = ident("py_const_" & $cnst)
  result = quote do:
    proc `procname`(): type(`cnst`) {.exportpy.} = `cnst`

macro exportdefaultconst(cnst: untyped): untyped =
  let procname = ident("py_const_" & $cnst)
  result = quote do:
    let py = pyBuiltinsModule()
    proc `procname`(): type(py.None) {.exportpy.} = py.None

template exportoptional(sym: untyped): untyped =
  when not declared(sym):
    exportdefaultconst(sym)
  else:
    exportconst(sym)

template export_plugin_metadata*(): untyped =
  exportoptional(ID)
  exportoptional(VERSION)
  exportoptional(INPUT)
  exportoptional(OUTPUT)
  exportoptional(PARAMETERS)
  exportoptional(METHOD)
  exportoptional(IMPLEMENTATION)
  exportoptional(REQ_SOFTWARE)
  exportoptional(REQ_HARDWARE)
  exportoptional(ADVICE)
