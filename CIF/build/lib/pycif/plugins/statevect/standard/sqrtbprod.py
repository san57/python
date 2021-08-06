import numpy as np


def sqrtbprod(statevect, chi, **kwargs):
    """Multiplies Chi by B**0.5.


    """

    # Initializes output vector
    xout = np.zeros(statevect.dim)

    # Loop over components of the control vector
    components = statevect.components
    for comp in components.attributes:
        component = getattr(components, comp)
        
        # Skip if component does not have parameters
        if not hasattr(component, "parameters"):
            continue
        
        for trcr in component.parameters.attributes:
            tracer = getattr(component.parameters, trcr)

            x_pointer = tracer.xpointer
            x_dim = tracer.dim
            chi_dim = tracer.chi_dim
            chi_pointer = tracer.chi_pointer
            ndates = tracer.ndates

            # Dealing with non-diagonal spatial B
            if tracer.hresol == "hpixels" and hasattr(tracer, "hcorrelations"):
                corr = tracer.hcorrelations
                sqrt_evalues = corr.sqrt_evalues
                evectors = corr.evectors

                # Re-stacking chi to period stacks
                chi_tempstacks = (
                    chi[chi_pointer: chi_pointer + chi_dim]
                    .transpose()
                    .reshape((tracer.chi_hresoldim, -1), order="F")
                )

                # Stack-multiplication of matrices
                # and flattening to control space
                chi_tmp = np.dot(
                    evectors, chi_tempstacks * sqrt_evalues[:, np.newaxis]
                ).flatten(order="F")

            else:
                chi_tmp = chi[chi_pointer: chi_pointer + chi_dim]

            # Deals with non-diagonal temporal correlations
            if hasattr(tracer, "tcorrelations"):
                corr = tracer.tcorrelations
                sqrt_evalues = corr.sqrt_evalues
                evectors = corr.evectors

                # Re-stacking chi
                chi_horizstacks = (
                    chi_tmp.transpose()
                    .reshape((-1, tracer.chi_tresoldim), order="F")
                    .T
                )

                # Stack-multiplication of matrices
                # and flattening to control space
                chi_tmp = np.dot(
                    evectors, chi_horizstacks * sqrt_evalues[:, np.newaxis]
                ).T.flatten(order="F")
            
            # TODO: deal with vertical resolution
            
            # Filling corresponding part in the control vector
            xout[x_pointer: x_pointer + x_dim] = chi_tmp

    return xout * statevect.std + statevect.xb


def sqrtbprod_ad(statevect, dx, **kwargs):
    # Initializes output vector
    chiout = np.zeros(statevect.chi_dim)

    # Loop over components of the control vector
    components = statevect.components
    for comp in components.attributes:
        component = getattr(components, comp)

        # Skip if component does not have parameters
        if not hasattr(component, "parameters"):
            continue
        
        for trcr in component.parameters.attributes:
            tracer = getattr(component.parameters, trcr)

            x_pointer = tracer.xpointer
            x_dim = tracer.dim
            chi_dim = tracer.chi_dim
            chi_pointer = tracer.chi_pointer
            chi_dim = tracer.chi_dim
            ndates = tracer.ndates

            # x * std
            xstd = (
                dx[x_pointer: x_pointer + x_dim]
                * statevect.std[x_pointer: x_pointer + x_dim]
            )

            # Dealing with non-diagonal temporal B
            if hasattr(tracer, "tcorrelations"):
                corr = tracer.tcorrelations
                sqrt_evalues = corr.sqrt_evalues
                evectors = corr.evectors

                # Re-stacking x to period stacks
                x_horizstacks = (
                    xstd.transpose()
                    .reshape((tracer.hresoldim, ndates), order="F")
                    .T
                )

                # Stack-multiplication of matrices
                # and flattening to control space
                xstd = np.dot(
                    evectors.T * sqrt_evalues[:, np.newaxis], x_horizstacks
                ).T.flatten(order="F")

            # Dealing with non-diagonal spatial B
            if tracer.hresol == "hpixels" and hasattr(tracer, "hcorrelations"):
                corr = tracer.hcorrelations
                sqrt_evalues = corr.sqrt_evalues
                evectors = corr.evectors

                # Re-stacking x to horizontal stacks
                x_tempstacks = xstd.transpose().reshape(
                    (tracer.hresoldim, tracer.chi_tresoldim), order="F"
                )

                # Stack-multiplication of matrices
                # and flattening to control space
                xstd = np.dot(
                    evectors.T * sqrt_evalues[:, np.newaxis], x_tempstacks
                ).flatten(order="F")

            # Filling Chi
            chiout[chi_pointer: chi_pointer + chi_dim] = xstd

    return chiout
