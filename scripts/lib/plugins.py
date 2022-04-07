COMPUTE_PLUGIN_INTERFACE = {}
COMPUTE_PLUGIN_INTERFACE["req_func"] = ["compute"]
COMPUTE_PLUGIN_INTERFACE["opt_func"] = ["initialize", "finalize"]
COMPUTE_PLUGIN_INTERFACE["req_const"] = ["ID", "VERSION", "INPUT", "OUTPUT"]
COMPUTE_PLUGIN_INTERFACE["opt_const"] = ["PARAMETERS", "METHOD",
                                         "IMPLEMENTATION", "ADVICE",
                                         "REQ_SOFTWARE", "REQ_HARDWARE"]

IDPROC_PLUGIN_INTERFACE = {}
IDPROC_PLUGIN_INTERFACE["req_func"] = ["compute_id"]

def plugin_metadata_str(plugin):
  metadata_keys = COMPLETE_PLUGIN_INTERFACE["req_const"] + \
                  COMPLETE_PLUGIN_INTERFACE["opt_const"]
  result = {k.lower(): getattr(plugin, k) \
              for k in metadata_keys if hasattr(plugin, k)}
  if "output" in result:
    result["output"] = ",".join(result["output"])
  if "parameters" in result:
    result["parameters"] = ";".join([",".join(e) for e in result["parameters"]])
  return result

