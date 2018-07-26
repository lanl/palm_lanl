def write_config(np_min, np_max, nx_min, nx_max, ny_min, ny_max, nz_min, nz_max, tpn, dnpexnpey, dnpexnpey_tolerance, dnxny, poisfft, switch,
                 temperton, mlt_grid, spectr, rslt_thrs, ld_thrs):
    #from configparser import ConfigParser  # windows
    from ConfigParser import ConfigParser  # Linux

    # if __name__ == "__main__":
    cfg = ConfigParser()

    import os.path
    if os.path.exists('.palm_gf_conf') is False:
        cfg.add_section("processor_topology")
        cfg.add_section("numerical_grid")
        cfg.add_section("method")

    #cfg.add_section("processor_topology")
    cfg.set("processor_topology", "np_min", np_min)
    cfg.set("processor_topology", "np_max", np_max)
    cfg.set("processor_topology", "tasks_per_node", tpn)
    cfg.set("processor_topology", "dnpexnpey", dnpexnpey)
    cfg.set("processor_topology", "dnpexnpey_tolerance", dnpexnpey_tolerance)
    cfg.set("processor_topology", "dnxny", dnxny)

    #cfgadd_section("numerical_grid").
    cfg.set("numerical_grid", "nx_min", nx_min)
    cfg.set("numerical_grid", "nx_max", nx_max)
    cfg.set("numerical_grid", "ny_min", ny_min)
    cfg.set("numerical_grid", "ny_max", ny_max)
    cfg.set("numerical_grid", "nz_min", nz_min)
    cfg.set("numerical_grid", "nz_max", nz_max)

    #cfg.add_section("method")
    cfg.set("method", "poisfft", poisfft)
    cfg.set("method", "switch", switch)
    cfg.set("method", "temperton", temperton)
    cfg.set("method", "mlt_grid", mlt_grid)
    cfg.set("method", "spectr", spectr)

    #cfg.add_section("settings")
    #cfg.set("settings", "path", path)
    #cfg.set("settings", "result_threshold", rslt_thrs)
    #cfg.set("settings", "load_threshold", ld_thrs)


    with open(".palm_gf_config", "a") as configfile:
        cfg.write(configfile)


def read_config():
    #from configparser import ConfigParser  # windows
    from ConfigParser import ConfigParser  # Linux


    cfg = ConfigParser()

    cfg.read(".palm_gf_config")
    np_min = cfg.get("processor_topology", "np_min")
    np_max = cfg.get("processor_topology", "np_max")
    tpn = cfg.get("processor_topology", "tasks_per_node")
    dnpexnpey = cfg.get("processor_topology", "dnpexnpey")
    dnpexnpey_tolerance = cfg.get("processor_topology", "dnpexnpey_tolerance")
    dnxny = cfg.get("processor_topology", "dnxny")

    nx_min = cfg.get("numerical_grid", "nx_min")
    nx_max = cfg.get("numerical_grid", "nx_max")
    ny_min = cfg.get("numerical_grid", "ny_min")
    ny_max = cfg.get("numerical_grid", "ny_max")
    nz_min = cfg.get("numerical_grid", "nz_min")
    nz_max = cfg.get("numerical_grid", "nz_max")

    poisfft = cfg.get("method", "poisfft")
    switch = cfg.get("method", "switch")
    temperton = cfg.get("method", "temperton")
    mlt_grid = cfg.get("method", "mlt_grid")
    spectr = cfg.get("method", "spectr")

    import ConfigParser as conf
    try:
        result_threshold = cfg.get("settings", "result_threshold")
        load_threshold = cfg.get("settings", "load_threshold")
        path = cfg.get("settings", 'path')

    except conf.NoSectionError:

        path = '/localdata/'
        result_threshold = 250000
        load_threshold = 100000

    return np_min, np_max, tpn, dnpexnpey, dnpexnpey_tolerance, dnxny, nx_min, nx_max, ny_min, ny_max, nz_min, nz_max, poisfft, switch, temperton, mlt_grid, spectr, result_threshold, load_threshold, path


def write_config_settings(path, rslt_thrs, ld_thrs):

    from ConfigParser import ConfigParser

    cfg = ConfigParser()

    cfg.add_section('settings')

    cfg.set('settings', 'path', path)
    cfg.set('settings', "result_threshold", rslt_thrs)
    cfg.set('settings', "load_threshold", ld_thrs)

    with open(".palm_gf_config", "a") as configfile:
        cfg.write(configfile)


def read_config_settings():

    from ConfigParser import ConfigParser
    import ConfigParser as con

    cfg = ConfigParser()
    cfg.read(".palm_gf_config")

    try:
        path = cfg.get("settings", "path")
        result_thrs = cfg.get("settings", "result_threshold")
        load_thrs = cfg.get("settings", "load_threshold")

    except con.NoSectionError:

        path = '/localdata/'
        result_thrs = 250000
        load_thrs = 100000

    except con.NoOptionError:

        path = '/localdata/'
        result_thrs = 250000
        load_thrs = 100000

    #print path, result_thrs, load_thrs

    return path, result_thrs, load_thrs


# ***********************************************


def closing_cleanup():
    from ConfigParser import ConfigParser  # Linux
    import ConfigParser as conf
    import os, shutil, time

    cfg = ConfigParser()

    try:

        cfg.read(".palm_gf_config")
        np_min = cfg.get("processor_topology", "np_min")
        np_max = cfg.get("processor_topology", "np_max")
        tpn = cfg.get("processor_topology", "tasks_per_node")
        dnpexnpey = cfg.get("processor_topology", "dnpexnpey")
        dnpexnpey_tolerance = cfg.get("processor_topology", "dnpexnpey_tolerance")
        dnxny = cfg.get("processor_topology", "dnxny")

        nx_min = cfg.get("numerical_grid", "nx_min")
        nx_max = cfg.get("numerical_grid", "nx_max")
        ny_min = cfg.get("numerical_grid", "ny_min")
        ny_max = cfg.get("numerical_grid", "ny_max")
        nz_min = cfg.get("numerical_grid", "nz_min")
        nz_max = cfg.get("numerical_grid", "nz_max")

        poisfft = cfg.get("method", "poisfft")
        switch = cfg.get("method", "switch")
        temperton = cfg.get("method", "temperton")
        mlt_grid = cfg.get("method", "mlt_grid")
        spectr = cfg.get("method", "spectr")

        var1_bool = True

    except conf.NoSectionError:
        np_min = 0
        np_max = 0
        tpn = 0
        dnpexnpey = 0
        dnpexnpey_tolerance = 0
        dnxny = 0

        nx_min = 0
        nx_max = 0
        ny_min = 0
        ny_max = 0
        nz_min = 0
        nz_max = 0

        poisfft = False
        switch = False
        temperton = False
        mlt_grid = False
        spectr = False
        var1_bool = False

    with open(".palm_gf_config", "w") as configfile:
        cfg.write(configfile)

    try:
        cfg.read(".palm_gf_config")
        rslt_thrs = cfg.get("settings", "result_threshold")
        ld_thrs = cfg.get("settings", "load_threshold")
        path = cfg.get("settings", 'path')

        var2_bool = True

    except conf.NoSectionError:

        path = '/localdata/'
        rslt_thrs = 250000
        ld_thrs = 100000

        var2_bool = False



    with open(".palm_gf_config", "w") as configfile:
        cfg.write(configfile)



    #os.remove('.palm_gf_config')



    try:
        cfg.set("processor_topology", "np_min", np_min)
        cfg.set("processor_topology", "np_max", np_max)
        cfg.set("processor_topology", "tasks_per_node", tpn)
        cfg.set("processor_topology", "dnpexnpey", dnpexnpey)
        cfg.set("processor_topology", "dnpexnpey_tolerance", dnpexnpey_tolerance)
        cfg.set("processor_topology", "dnxny", dnxny)


        cfg.set("numerical_grid", "nx_min", nx_min)
        cfg.set("numerical_grid", "nx_max", nx_max)
        cfg.set("numerical_grid", "ny_min", ny_min)
        cfg.set("numerical_grid", "ny_max", ny_max)
        cfg.set("numerical_grid", "nz_min", nz_min)
        cfg.set("numerical_grid", "nz_max", nz_max)


        cfg.set("method", "poisfft", poisfft)
        cfg.set("method", "switch", switch)
        cfg.set("method", "temperton", temperton)
        cfg.set("method", "mlt_grid", mlt_grid)
        cfg.set("method", "spectr", spectr)

    except conf.NoSectionError:
        cfg.add_section("processor_topology")
        cfg.set("processor_topology", "np_min", np_min)
        cfg.set("processor_topology", "np_max", np_max)
        cfg.set("processor_topology", "tasks_per_node", tpn)
        cfg.set("processor_topology", "dnpexnpey", dnpexnpey)
        cfg.set("processor_topology", "dnpexnpey_tolerance", dnpexnpey_tolerance)
        cfg.set("processor_topology", "dnxny", dnxny)

        cfg.add_section("numerical_grid")
        cfg.set("numerical_grid", "nx_min", nx_min)
        cfg.set("numerical_grid", "nx_max", nx_max)
        cfg.set("numerical_grid", "ny_min", ny_min)
        cfg.set("numerical_grid", "ny_max", ny_max)
        cfg.set("numerical_grid", "nz_min", nz_min)
        cfg.set("numerical_grid", "nz_max", nz_max)

        cfg.add_section("method")
        cfg.set("method", "poisfft", poisfft)
        cfg.set("method", "switch", switch)
        cfg.set("method", "temperton", temperton)
        cfg.set("method", "mlt_grid", mlt_grid)
        cfg.set("method", "spectr", spectr)


    try:
        cfg.set("settings", "path", path)
        cfg.set("settings", "result_threshold", rslt_thrs)
        cfg.set("settings", "load_threshold", ld_thrs)

    except conf.NoSectionError:
        cfg.add_section("settings")
        cfg.set("settings", "path", path)
        cfg.set("settings", "result_threshold", rslt_thrs)
        cfg.set("settings", "load_threshold", ld_thrs)






    with open(".palm_gf_config", "w") as configfile:
        cfg.write(configfile)

