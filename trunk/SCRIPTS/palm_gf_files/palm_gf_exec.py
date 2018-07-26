def grid_executer():
    import palm_gf_conf as config_wr
    import sqlite3
    conn = sqlite3.connect(".palm_gf_data.db")
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS grid_current")
    c.execute("DROP TABLE IF EXISTS grid_limits")
    c.execute("CREATE TABLE IF NOT EXISTS grid_current(nx INT, ny INT, nz INT, npex INT, npey INT, npxnpy FLOAT, np INT, ngpts INT)")
    c.execute("CREATE TABLE IF NOT EXISTS grid_limits(nx INT, ny INT, nz INT, npex INT, npey INT, npxnpy FLOAT, np INT, ngpts INT)")
    conn.commit()
    main_bool = True
    parameters = config_wr.read_config()
    print(parameters)
    min_procs = int(parameters[0])
    max_procs = int(parameters[1])
    tpn = int(parameters[2])
    dnpexnpey = float(parameters[3])
    dnxny = float(parameters[4])
    nx_min = int(parameters[5])
    nx_max = int(parameters[6])
    ny_min = int(parameters[7])
    ny_max = int(parameters[8])
    nz_min = int(parameters[9])
    nz_max = int(parameters[10])
    poisfft = parameters[11]
    switch = parameters[12]
    tempterton = parameters[13]
    mlt_grid = parameters[14]
    spectr = parameters[15]

    if poisfft == str(True):
        poisfft = True
    else:
        poisfft = False

    if switch == str(True):
        switch = True
    else:
        switch = False

    if tempterton == str(True):
        tempterton = True
    else:
        tempterton = False

    if mlt_grid == str(True):
        mlt_grid = True
    else:
        mlt_grid = False

    if spectr == str(True):
        spectr = True
    else:
        spectr = False

    print(spectr, type(spectr))
    print(poisfft, switch, tempterton, mlt_grid, spectr)

    np_used = min_procs
    counter = 0

    nx = nx_min
    ny = ny_min
    nz = nz_min

    from math import floor

    def factors(n):
        result = []
        for i in range(2, n + 1):  # test all integers between 2 and n
            s = 0
            while n / i == floor(n / float(i)):  # is n/i an integer?
                n = n / float(i)
                s += 1
            if s > 0:
                for k in range(s):
                    result.append(i)  # i is a pf s times
                    if n == 1:
                        return result

    while np_used <= max_procs:
        a = 1

        while a <= np_used:
            prcs_var = np_used % a
            if prcs_var != 0:
                a += 1
            elif prcs_var == 0:
                npex = a
                npey = int(np_used / npex)

                if tpn != 0:
                    if np_used % tpn != 0:
                        a += 1
                        continue

                if dnpexnpey != 0 and npex / npey != dnpexnpey:
                    a += 1
                    continue

                while nx <= nx_max:
                    if (nx + 1) % npex != 0:
                        nx += 1
                        continue

                    if mlt_grid is True and (nx + 1) % 2 != 0:
                        nx += 1
                        continue

                    if switch is True and (nx + 1) % npey != 0:
                        nx += 1
                        continue
                    if npex > nx:
                        nx += 1
                        continue

                    while ny <= ny_max:

                        if dnxny != 0 and float(nx) / float(ny) != float(dnxny):
                            ny += 1
                            continue
                        if (ny + 1) % npey != 0:
                            ny += 1
                            continue

                        if mlt_grid is True and ny % 2 != 0:
                            ny += 1
                            continue

                        if (ny + 1) % npex != 0 and switch is True:
                            ny += 1
                            continue
                        if npey > ny:
                            ny += 1
                            continue

                        while nz <= nz_max:

                            if mlt_grid is True and nz % 2 != 0:
                                nz += 1
                                continue

                            if poisfft is True and nz % npex != 0:
                                nz += 1
                                continue

                            if spectr is True and nz % npey != 0:
                                nz += 1
                                continue

                            if tempterton is True and nx > 1 and ny > 1:  # and nz < 1:

                                nx_list = factors(nx + 1)
                                i = 0
                                nx_var = nx_list[i]

                                while i < len(nx_list):
                                    if nx_var != 2 or nx_var != 3 or nx_var != 5:

                                        i += 1
                                        continue

                                    i += 1
                                ny_list = factors(ny + 1)
                                i = 0
                                ny_var = ny_list[i]
                                while i < len(ny_list):
                                    if ny_var != 2 or ny_var != 3 or ny_var != 5:

                                        i += 1
                                        continue
                                    i += 1

                            counter += 1

                            npxnpy = format(float(npex)/float(npey), '.2f')
                            c.execute("""INSERT OR REPLACE INTO grid_current(nx, ny, nz, npex, npey, npxnpy, np, ngpts) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                            (nx, ny, nz, npex, npey, npxnpy, (npex * npey), (nx*ny*nz)))

                            nz += 11
                        nz = nz_min
                        ny += 1
                    ny = ny_min
                    nx += 1
                nx = nx_min
                a += 1
                #  a += 1
        np_used += 1

        conn.commit()



    conn.commit()
    c.close()
    conn.close()

    #********************************

    conn = sqlite3.connect(".palm_gf_data.db")
    c = conn.cursor()
    try:
        c.execute("SELECT nx FROM grid_current ORDER BY nx DESC LIMIT 1")
        mx_nx = c.fetchone()[0]
        #print(mx_nx)
        c.execute("SELECT nx FROM grid_current ORDER BY nx  LIMIT 1")
        mn_nx = c.fetchone()[0]
        #print(mn_nx)
        c.execute("SELECT ny FROM grid_current ORDER BY ny DESC LIMIT 1")
        mx_ny = c.fetchone()[0]
        #print(mx_ny)
        c.execute("SELECT ny FROM grid_current ORDER BY ny  LIMIT 1")
        mn_ny = c.fetchone()[0]
        #print(mn_ny)
        c.execute("SELECT nz FROM grid_current ORDER BY nz DESC LIMIT 1")
        mx_nz = c.fetchone()[0]
        #print(mx_nz)
        c.execute("SELECT nz FROM grid_current ORDER BY nz  LIMIT 1")
        mn_nz = c.fetchone()[0]
        #print(mn_nz)
        c.execute("SELECT npex FROM grid_current ORDER BY npex DESC LIMIT 1")
        mx_npex = c.fetchone()[0]
        #print(mx_npex)
        c.execute("SELECT npex FROM grid_current ORDER BY npex  LIMIT 1")
        mn_npex = c.fetchone()[0]
        #print(mn_npex)
        c.execute("SELECT npey FROM grid_current ORDER BY npey DESC LIMIT 1")
        mx_npey = c.fetchone()[0]
        #print(mx_npey)
        c.execute("SELECT npey FROM grid_current ORDER BY npey  LIMIT 1")
        mn_npey = c.fetchone()[0]
        #print(mn_npey)
        c.execute("SELECT npxnpy FROM grid_current ORDER BY npxnpy DESC LIMIT 1")
        mx_npxnpy = c.fetchone()[0]
        #print(mx_npxnpy)
        c.execute("SELECT npxnpy FROM grid_current ORDER BY npxnpy  LIMIT 1")
        mn_npxnpy = c.fetchone()[0]
        #print(mn_npxnpy)
        c.execute("SELECT np FROM grid_current ORDER BY np DESC LIMIT 1")
        mx_np = c.fetchone()[0]
        #print(mx_np)
        c.execute("SELECT np FROM grid_current ORDER BY np LIMIT 1")
        mn_np = c.fetchone()[0]
        #print(mn_np)
        c.execute("SELECT ngpts FROM grid_current ORDER BY ngpts DESC LIMIT 1")
        mx_ngpts = c.fetchone()[0]
        #print(mx_ngpts)
        c.execute("SELECT ngpts FROM grid_current ORDER BY ngpts LIMIT 1")
        mn_ngpts = c.fetchone()[0]
        #print(mn_ngpts)

        conn.commit()
        c.execute(
            """INSERT OR REPLACE INTO grid_limits(nx, ny, nz, npex, npey, npxnpy, np, ngpts) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (mn_nx, mn_ny, mn_nz, mn_npex, mn_npey, mn_npxnpy, mn_np, mn_ngpts))
        #conn.commit()
        #print("test", mn_nx)
        c.execute(
            """INSERT OR REPLACE INTO grid_limits(nx, ny, nz, npex, npey, npxnpy, np, ngpts) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (mx_nx, mx_ny, mx_nz, mx_npex, mx_npey, mx_npxnpy, mx_np, mx_ngpts))
        conn.commit()

        c.close()
        conn.close()
    except TypeError:

        checkfile = open(".palm_gf_tmp", "w")
        if counter != 0:
            checkfile.write("Gridfinder found " + str(counter) + " results.\n1")
        else:
            checkfile.write("Check input, no Results found.\n0")
        checkfile.close()



    checkfile = open(".palm_gf_tmp", "w")
    if counter != 0:
        checkfile.write("Gridfinder found " + str(counter) + " results.\n1")
    else:
        checkfile.write("Check input, no Results found.\n0")
    checkfile.close()
    return main_bool



def grid_maxmin():

    import palm_gf_conf
    import sqlite3

    conn = sqlite3.connect("grid_current.db")
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS grid_maxmin")
    c.execute("CREATE TABLE IF NOT EXISTS grid_maxmin(np_min INT,np_max INT, npex_min INT,npex_max INT, npey_min INT,npey_max INT, nx_min INT,nx_max INT,"
              " ny_min INT,ny_max INT, nz_min INT, nz_max INT)")
    parm = config_wr.read_config()
    min_procs = int(parm[0])
    max_procs = int(parameters[1])
    tpn = int(parameters[2])
    dnpexnpey = float(parameters[3])
    dnxny = float(parameters[4])
    nx_min = int(parameters[5])
    nx_max = int(parameters[6])
    ny_min = int(parameters[7])
    ny_max = int(parameters[8])
    nz_min = int(parameters[9])
    nz_max = int(parameters[10])





    c.execute(
        """INSERT OR REPLACE INTO grid_maxmin(np_min, np_max, npex_min, npex_max, npey_min, npey_max, nx_min, nx_max,
         ny_min, ny_max, nz_min, nz_max) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

    conn.commit()
    c.close()
    conn.close()