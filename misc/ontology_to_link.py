#!/usr/bin/env python3
import yaml

# do this once only by setting up the application:
with open("ontology_links.yaml") as f:
  ontology_links_table = yaml.safe_load(f)

#
# Returns None if no hyperlink can be computed, in which case
# no hyperlink shall be visualized
#
def get_ontology_hyperlink(ontology_ref, ontology_links_table):
  ontology_pfx, term_id = ontology_ref.split(":", 1)
  if ontology_pfx in ontology_links_table:
    return ontology_links_table[ontology_pfx].format(ID=term_id)
  else:
    return None

def handle_ontology_xref(ontology_xref, ontology_links_table):
  if not ontology_xref:
    # do not show anything if there is no content
    return
  else:
    label = "Ontology term: "
    hyperlink = get_ontology_hyperlink(ontology_xref, ontology_links_table)
    if hyperlink is None:
      print(label + ontology_xref)
    else:
      # in the web app use a hyperlink under ontology_xref instead, here --> is
      # just to demonstrate:
      print(label + ontology_xref + " --> " + hyperlink)

def handle_related_ontology_terms(related_ontology_terms, ontology_links_table):
  if not related_ontology_terms:
    # do not show anything if there is no content
    return
  else:
    label = "Related ontology terms:"
    content = [label]
    for line in related_ontology_terms.split("\n"):
      pfx, ontology_xref = line.split(": ", 1)
      hyperlink = get_ontology_hyperlink(ontology_xref, ontology_links_table)
      if hyperlink is None:
        content.append(f"- {ontology_xref}")
      else:
        # in the web app use a hyperlink under ontology_xref instead, here --> is
        # just to demonstrate:
        content.append(f"- {ontology_xref} --> {hyperlink}")
    print("\n".join(content))

# this simulates the results of a query

for ontology_xref in ["SO:0001459", "EDAM_data:129", ""]:
  handle_ontology_xref(ontology_xref, ontology_links_table)
for related_ontology_terms in ["", "genome: SO:0001026", "genome: SO:0001026\nchromosome: SO:0000340\nplasmid: SO:0000155\nreplicon: SO:0001235"]:
  handle_related_ontology_terms(related_ontology_terms, ontology_links_table)
