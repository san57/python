from pycif.utils.path import init_dir
from .fetch import default_fetch


def init_components(plugin):
    if hasattr(plugin, "components"):
        components = plugin.components
        for comp in components.attributes:
            component = getattr(components, comp)

            # Fetch parameters
            # If no parameters, handle the component as a whole
            if not hasattr(component, "parameters"):
                params = component
                parameters = [""]
            else:
                params = component.parameters
                parameters = params.attributes[:] + [""]

            # Loop over parameters to fetch information
            for trcr in parameters:
                tracer = getattr(params, trcr, component)

                # By default is not in the target vector
                tracer.iscontrol = False

                # Fetch reference directory and file format
                trac_dir = getattr(
                    tracer,
                    "dir",
                    getattr(component, "dir", getattr(components, "dir", "")),
                )
                trac_file = getattr(
                    tracer, "file", getattr(component, "file", "")
                )
                tracer.varname = getattr(tracer, "varname", "")

                # Initializes target directory and pass info to tracer
                target_dir = "{}/statevect/{}/{}".format(
                    plugin.workdir, comp, trcr
                )
                init_dir(target_dir)

                tracer.dir = target_dir
                tracer.file = trac_file

                # Forces the tracer to have an empty read function
                if not hasattr(tracer, "read"):
                    tracer = tracer.get_subclass("fields")(plg_orig=tracer)
                    if trcr != "":
                        setattr(component.parameters, trcr, tracer)
                    else:
                        setattr(plugin.components, comp, tracer)

                # Gets read/fetch from model if not already defined
                # Passing the tracer to the read function for attribute access
                mod = plugin.model
                backups = list(mod.backup_comps.values())
                cmp = (
                    comp
                    if hasattr(mod, comp)
                    else (list(mod.backup_comps.keys())[backups.index(comp)]
                          if comp in backups else None)
                    )
                if cmp is None:
                    raise Exception("{} in your Yaml is not recognized as"
                                    " a valid input for the model".format(comp))
                
                mod_comp = getattr(plugin.model, cmp)

                tracer.read = getattr(tracer, "read", mod_comp.read)
                tracer.fetch = getattr(
                    tracer, "fetch", getattr(mod_comp, "fetch", default_fetch)
                )

                for attr in ["file", "dir"]:
                    if getattr(tracer, attr, "") == "":
                        setattr(tracer, attr, getattr(mod_comp, attr, ""))

                # Fetch files and dates
                list_files, list_dates = tracer.fetch(
                    trac_dir,
                    trac_file,
                    plugin.model.input_dates,
                    target_dir,
                    component=component,
                    tracer=tracer,
                )
                tracer.input_dates = list_dates
                tracer.input_files = list_files

                # Saving tracer_dir into component dir if not already available
                if not hasattr(component, "input_dates"):
                    component.input_dates = list_dates
                if not hasattr(component, "input_files"):
                    component.input_files = list_files

                # Get the domain and
                # change it to the domain side if lateral conditions
                if hasattr(tracer, "domain"):
                    continue

                if hasattr(tracer, "get_domain"):
                    tracer.domain = tracer.get_domain(
                        trac_dir,
                        trac_file,
                        plugin.model.input_dates,
                        target_dir,
                        tracer=tracer,
                    )

                else:
                    tracer.domain = plugin.domain
