# Interface of a plugin

# Compute function

The plugin shall export the `compute(filename, **kwargs)` function:
- `filename`: name of the file with the input data
- `kwargs`: additional optional named parameters (described in the PARAMETERS
            constant, see "Metadata constants section" below)

The return value of the method is a 2-tuple `(results, logs)`
where:
 - `results`: list of values, in the order specified by
              the `OUTPUT` constant (see below),
 - `logs`: a possibly empty strings list, containing messages to be displayed
           to the user or stored in a log file;
           suggested format: "{key}\t{message}",
           where key identifies the type of log message.

## Metadata constants

The plugin communicates its purpose, version and interface by defining
the following constants.

Mandatory (their combination must be unique among all plugins):
 - `ID`:      name of the plugin              [str, max len 256]
 - `VERSION`: a version number or name        [str, max len 64]

Mandatory:
 - `INPUT`:   short description of the input (type of data, format, source)
              [str, max len 512]
 - `OUTPUT`:  list of strings (attribute ids), describing the order or the
              attributes in the output; each of them must be defined in
              `attribute_definitions.yaml`.

Required if compute() accepts named parameters:
 - `PARAMETERS`: optional named parameters of the compute function;
                 list of 4-tuples of strings:
                 (name, datatype, default_value, documentation)

Optional: (string, max len 4096)
 - `METHOD`:         how the results are computed, conceptually
 - `IMPLEMENTATION`: how the method is implemented, technically
 - `REQ_SOFTWARE`:   required tools or libraries
 - `REQ_HARDWARE`:   required hardware resources (memory, GPUs...)
 - `ADVICE`:         when should this method used instead of others

## Nim plugins

Plugins can be implemented in Nim using the `nimpy` library.

Besides `nimpy`, the plugin shall import the `plugins_helpers` module
and call the template `export_template_metadata()`, after the definition
of the metadata constants.

### Compute function in Nim

The compute function has the following signature, if there are no
optional parameters:

```
proc compute(filename: string):
             tuple[results: seq[string]], logs: seq[string]] {.exportpy.}
```

Any accepted optional parameter must be defined in the signature.
E.g. given two optional named parameters p1 of type t1 and default value d1
and p2 of type t2 and default value d2, the signature becomes:
```
proc compute(filename: string, p1: t1 = d1, p2: t2 = d2):
             tuple[results: seq[string]], logs: seq[string]] {.exportpy.}
```

### Internals of constant export

The `nimpy` library does not yet allow to export constants to Python.

The `export_template_metadata` template exports all plugin-related constants as
procs with a prefix `py_const_` and no arguments, returning the value of the
constant (e.g. `py_const_PARAMETERS()`) or a Python None value if an optional
constant was not defined.

The python library `mod.py` under `lib` stores the values returned by
`py_const_` in module variables without the prefix (e.g. `PARAMETERS`) and
removes the prefixed procs from the module `__dict__`.

