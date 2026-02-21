import functools
import os
import shutil

from .log import logger


def latex_plot_format(func):
    """
    Wrap the plotting function to set rcParams dependent on environment variables

    The rcparams can be set directly from the env. variable `BILBY_STYLE` to
    point to a matplotlib style file. Or, if `BILBY_STYLE=default` (any case) a
    default setup is used, this is enabled by default. To not use any rcParams,
    set `BILBY_STYLE=none`. Occasionally, issues arrise with the latex
    `mathdefault` command. A fix is to define this command in the rcParams. An
    env. variable `BILBY_MATHDEFAULT` can be used to turn this fix on/off.
    Setting `BILBY_MATHDEFAULT=1` will enable the fix, all other choices
    (including undefined) will disable it. Additionally, the BILBY_STYLE and
    BILBY_MATHDEFAULT arguments can be passed into any
    latex_plot_format-wrapped plotting function and will be set directly.

    """
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib import rcParams

        if "BILBY_STYLE" in kwargs:
            bilby_style = kwargs.pop("BILBY_STYLE")
        else:
            bilby_style = os.environ.get("BILBY_STYLE", "default")

        if "BILBY_MATHDEFAULT" in kwargs:
            bilby_mathdefault = kwargs.pop("BILBY_MATHDEFAULT")
        else:
            bilby_mathdefault = int(os.environ.get("BILBY_MATHDEFAULT", "0"))

        if bilby_mathdefault == 1:
            logger.debug("Setting mathdefault in the rcParams")
            rcParams['text.latex.preamble'] = r'\providecommand{\mathdefault}[1][]{}'

        logger.debug("Using BILBY_STYLE={}".format(bilby_style))
        if bilby_style.lower() == "none":
            return func(*args, **kwargs)
        elif os.path.isfile(bilby_style):
            plt.style.use(bilby_style)
            return func(*args, **kwargs)
        elif bilby_style in plt.style.available:
            plt.style.use(bilby_style)
            return func(*args, **kwargs)
        elif bilby_style.lower() == "default":
            _old_tex = rcParams["text.usetex"]
            _old_serif = rcParams["font.serif"]
            _old_family = rcParams["font.family"]
            if shutil.which("latex"):
                rcParams["text.usetex"] = True
            else:
                rcParams["text.usetex"] = False
            rcParams["font.serif"] = "Computer Modern Roman"
            rcParams["font.family"] = "serif"
            rcParams["text.usetex"] = _old_tex
            rcParams["font.serif"] = _old_serif
            rcParams["font.family"] = _old_family
            return func(*args, **kwargs)
        else:
            logger.debug(
                "Environment variable BILBY_STYLE={} not used"
                .format(bilby_style)
            )
            return func(*args, **kwargs)
    return wrapper_decorator

def build_truth_dict(result):
    """
    Construct a dictionary of truth values for plotting or diagnostic purposes.

    The function extracts injection parameters from a result object and
    populates a dictionary of truth values corresponding to the search
    parameters. If component mass parameters are available, derived compact
    binary parameters such as chirp mass, mass ratio, and symmetric mass
    ratio are computed.

    Supported derived parameters
    ----------------------------
    - chirp_mass
        ((m1 * m2)^(3/5)) / (m1 + m2)^(1/5)

    - mass_ratio (or q)
        min(m1, m2) / max(m1, m2)

    - eta (symmetric mass ratio)
        (m1 * m2) / (m1 + m2)^2

    Parameters
    ----------
    result : object
        Result object containing inference outputs.

        Expected attributes:
        - injection_parameters : dict-like
            Dictionary of injected physical parameters.
        - search_parameter_keys : iterable
            List of parameters used in sampling/search.

    Returns
    -------
    truth_dict : dict
        Dictionary mapping parameter names to truth values.

        Contains:
        - Directly injected parameter values when available.
        - Derived compact-binary parameters when component masses
          are present.
        - NaN is used when a parameter cannot be determined.

    Notes
    -----
    - Mass-derived parameters are computed only if both component masses
      (mass_1 and mass_2) are available.
    - The function is designed to be compatible with posterior result
      objects stored in JSON or HDF5 formats, provided they are loaded
      into the appropriate Python object representation.

    Examples
    --------
    >>> truth_dict = build_truth_dict(result)
    >>> print(truth_dict['chirp_mass'])

    """

    inj = result.injection_parameters
    search = result.search_parameter_keys

    truth_dict = {}

    # Precompute masses once (if available)
    m1 = inj.get("mass_1")
    m2 = inj.get("mass_2")

    for p in search:

        # Direct injection parameter
        if p in inj:
            truth_dict[p] = inj[p]

        # Derived: chirp mass
        elif p == "chirp_mass" and m1 is not None and m2 is not None:
            truth_dict[p] = ((m1 * m2) ** (3/5)) / ((m1 + m2) ** (1/5))

        # Derived: mass ratio q (<= 1 convention)
        elif p in ["mass_ratio", "q"] and m1 is not None and m2 is not None:
            truth_dict[p] = min(m1, m2) / max(m1, m2)

        # Derived: symmetric mass ratio eta
        elif p in ["symmetric_mass_ratio", "eta"] and m1 is not None and m2 is not None:
            truth_dict[p] = (m1 * m2) / (m1 + m2)**2

    return truth_dict
