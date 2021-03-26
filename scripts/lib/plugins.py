PLUGIN_METADATA_KEYS = ["ID", "VERSION", "INPUT", "OUTPUT", "PARAMETER",
                        "METHOD", "IMPLEMENTATION", "REQ_SOFTWARE",
                        "REQ_HARDWARE", "ADVICE"]

def metadata(plugin):
  result = {k.lower(): getattr(plugin, k) for \
      k in PLUGIN_METADATA_KEYS if hasattr(plugin, k)}
  if "output" in result:
    result["output"] = ",".join(result["output"])
  return result
