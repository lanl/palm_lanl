import sys
import os
import sqlite3
from PyQt4 import QtCore, QtGui, uic
import subprocess as sub
import palm_gf_conf as configwr
   
out = sub.check_output("echo $PALM_BIN", shell=True, stderr=sub.STDOUT)
palm_bin = out.rstrip()

qtCreatorFile = palm_bin + '/palm_gf_files/palm_gf_table.ui'

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp3(QtGui.QMainWindow, Ui_MainWindow):

    class MyTableWidgetItem(QtGui.QTableWidgetItem):
        def __init__(self, text, sortKey):
            QtGui.QTableWidgetItem.__init__(self, text, QtGui.QTableWidgetItem.UserType)
            self.sortKey = sortKey

        def __lt__(self, other):
            return self.sortKey < other.sortKey

    def __init__(self):   #     def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)

        framegm = self.frameGeometry()
        screen = QtGui.QApplication.desktop().screenNumber(QtGui.QApplication.desktop().cursor().pos())
        centerpoint = QtGui.QApplication.desktop().screenGeometry(screen).center()
        framegm.moveCenter(centerpoint)

        #centerpoint = str(centerpoint)
        #xcenter = centerpoint.split('(')[1].split(',')[0]
        #ycenter = centerpoint.split('(')[1].split(',')[1].split(')')[0]
        ##print xcenter, ycenter
        #centerpoint = QtCore.QPoint(int(xcenter) + 418, int(ycenter))
        #framegm.moveCenter(centerpoint)
        self.move(framegm.topLeft())



        self.setupUi(self)
        self.pushButton_3.clicked.connect(self.load_result)
        #self.nx_min.valueChanged.connect(self.check)
        self.load_trigger()
        self.pushButton.clicked.connect(self.filter_results)
        #self.Sortnow_button.clicked.connect(self.sort_order)
        self.pushButton_2.clicked.connect(self.load_trigger)
        self.tableWidget.horizontalHeader().setClickable(True)


        self.nx_table.clicked.connect(lambda: self.sort_table(str("nx")))
        self.ny_table.clicked.connect(lambda: self.sort_table(str("ny")))
        self.nz_table.clicked.connect(lambda: self.sort_table(str("nz")))
        self.npex_table.clicked.connect(lambda: self.sort_table(str("npex")))
        self.npey_table.clicked.connect(lambda: self.sort_table(str("npey")))
        self.npexnpey_table.clicked.connect(lambda: self.sort_table(str("npexnpey")))
        self.np_table.clicked.connect(lambda: self.sort_table(str("np")))
        self.ngpts_table.clicked.connect(lambda: self.sort_table(str("ngpts")))
        self.nxpex_table.clicked.connect(lambda: self.sort_table(str("nxnpex")))
        self.nypey_table.clicked.connect(lambda: self.sort_table(str("nynpey")))

        self.instant()

        self.save_to_file_button.clicked.connect(self.get_path)

    def instant(self):
        checkfile = open(".palm_gf_tmp", "r")
        # self.result_label.setText(str(checkfile.readline()))
        res_text = str(checkfile.readline())  # XXX
        result_nr = res_text.split(' ')[2]
        checkfile.close()
        if int(result_nr) < 100000:
            self.load_result()



    def load_trigger(self):
        pathx = configwr.read_config()
        pathx = pathx[19]


        dtb = str('.palm_gf_data.db')
        #con = sqlite3.connect("/localdata/.palm_gf_data.db")

        pathx = pathx + '/.palm_gf_data.db'
        con = sqlite3.connect(pathx)
        c = con.cursor()
        c.execute("SELECT * FROM " + 'grid_limits')
        mini = c.fetchone()
        max = c.fetchone()
        self.nx_min.setValue(mini[0])
        self.nx_max.setValue(max[0])

        self.ny_min.setValue(mini[1])
        self.ny_max.setValue(max[1])
        self.nz_min.setValue(mini[2])
        self.nz_max.setValue(max[2])
        self.npex_min.setValue(mini[3])
        self.npex_max.setValue(max[3])
        self.npey_min.setValue(mini[4])
        self.npey_max.setValue(max[4])
        self.npxnpy_min.setValue(mini[5])
        self.npxnpy_max.setValue(max[5])
        self.np_min.setValue(mini[6])
        self.np_max.setValue(max[6])
        self.ngpts_min.setValue(mini[7])
        self.ngpts_max.setValue(max[7])
        self.nxpex_min.setValue(mini[8])
        self.nxpex_max.setValue(max[8])
        self.nypey_min.setValue(mini[9])
        self.nypey_max.setValue(max[9])

        self.nx_min.setMinimum(mini[0])
        self.nx_max.setMaximum(max[0])
        self.ny_min.setMinimum(mini[1])
        self.ny_max.setMaximum(max[1])
        self.nz_min.setMinimum(mini[2])
        self.nz_max.setMaximum(max[2])
        self.npex_min.setMinimum(mini[3])
        self.npex_max.setMaximum(max[3])
        self.npey_min.setMinimum(mini[4])
        self.npey_max.setMaximum(max[4])
        self.npxnpy_min.setMinimum(mini[5])
        self.npxnpy_max.setMaximum(max[5])
        self.np_min.setMinimum(mini[6])
        self.np_max.setMaximum(max[6])

        self.ngpts_min.setMinimum(mini[7])
        self.ngpts_max.setMaximum(max[7])
        self.ngpts_min.setMaximum(max[7])
        self.ngpts_max.setMinimum(mini[7])

        self.nxpex_min.setMinimum(mini[8])
        self.nxpex_max.setMaximum(max[8])
        self.nxpex_min.setMaximum(max[8])
        self.nxpex_max.setMinimum(mini[8])

        self.nypey_min.setMinimum(mini[9])
        self.nypey_max.setMaximum(max[9])
        self.nypey_min.setMaximum(max[9])
        self.nypey_max.setMinimum(mini[9])





        con.commit()
        c.close()
        con.close()




    def check(self):
        pathx = configwr.read_config()
        pathx = pathx[19]

        dtb = str('.palm_gf_data.db')
        #con = sqlite3.connect("/localdata/.palm_gf_data.db")
        con = sqlite3.connect(pathx + '/.palm_gf_data.db')
        c = con.cursor()
        c.execute("SELECT * FROM " + 'grid_limits')
        mini = c.fetchone()
        max = c.fetchone()
        if self.nx_min.value() < mini[0]:

            self.nx_min.setValue(mini[0])


    def process1(self):
        self.calc_label.setText('loading...')
        QtGui.QApplication.processEvents()



    def load_result(self):
        #print("LOADED!!!")
        import decimal

        pathx = configwr.read_config()
        pathx = pathx[19]

        self.setEnabled(False)
        QtGui.QApplication.processEvents()
        self.load_trigger()
        self.process1()
        database = str('.palm_gf_data.db')
        conn = sqlite3.connect(pathx + '/.palm_gf_data.db')
        c = conn.cursor()
        c.execute("SELECT * FROM " + 'grid_current')
        results = c.fetchall()

        self.tableWidget.setRowCount(len(results))

        i = 0
        j = 0
        k = 0

        while i < len(results):
            line = results[i]
            while j < 10:
                var = line[j]

                if j == 7:
                    self.tableWidget.setItem(i, j, self.MyTableWidgetItem(str("%.1e" % var), j + i))
                    #print("%.2e" % int(var), "%.4e" % int(1782))

                else:
                    self.tableWidget.setItem(i, j, self.MyTableWidgetItem(str(var), j+i))
                #item = self.MyTableWidgetItem(str(var), k)
                #self.tableWidget.setItem(i, j, item)

                j += 1
                #if j == 3:
                    #print(k)
                #k += 1
            #k -= 7
            k += 1
            j = 0
            i += 1

        c.close()
        conn.close()
        self.calc_label.setText('loading completed')
        self.label_11.setText(str(len(results)) + ' results ')
        self.setEnabled(True)
        QtGui.QApplication.processEvents()

    def filter_results(self):

        pathx = configwr.read_config()
        pathx = pathx[19]

        self.setEnabled(False)
        self.calc_label.setText('calculating...')
        QtGui.QApplication.processEvents()
        database = str('.palm_gf_data.db')
        conn = sqlite3.connect(pathx + '/.palm_gf_data.db')
        c = conn.cursor()
        c.execute("SELECT * FROM " + "grid_current")
        results = c.fetchall()
        #print(results)

        self.tableWidget.setRowCount(len(results))

        i = 0
        j = 0
        row_cnt = -1
        while i < len(results):
            line = results[i]

            if line[0] <= self.nx_max.value():

                if line[0] >= self.nx_min.value():

                    if line[1] <= self.ny_max.value():

                        if line[1] >= self.ny_min.value():

                            if line[2] <= self.nz_max.value():

                                if line[2] >= self.nz_min.value():

                                    if line[3] <= self.npex_max.value():

                                        if line[3] >= self.npex_min.value():

                                            if line[4] <= self.npey_max.value():

                                                if line[4] >= self.npey_min.value():

                                                    if line[5] <= self.npxnpy_max.value():

                                                        if line[5] >= self.npxnpy_min.value():

                                                            if line[6] <= self.np_max.value():

                                                                if line[6] >= self.np_min.value():

                                                                    if line[7] <= self.ngpts_max.value():

                                                                        if line[7] >= self.ngpts_min.value():

                                                                            if line[8] <= self.nxpex_max.value():

                                                                                if line[8] >= self.nxpex_min.value():

                                                                                    if line[9] <= self.nypey_max.value():

                                                                                        if line[9] >= self.nypey_min.value():

                                                                                            row_cnt += 1
                                                                                            while j < 10:
                                                                                                var = line[j]

                                                                                                if j == 7:
                                                                                                    self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(str("%.1e" % var), i))

                                                                                                else:
                                                                                                    self.tableWidget.setItem(row_cnt, j, self.MyTableWidgetItem(str(var), i))

                                                                                                j += 1

            j = 0
            i += 1

        c.close()
        conn.close()
        self.tableWidget.setRowCount(row_cnt + 1)
        self.setEnabled(True)
        self.calc_label.setText('calculation completed')
        self.label_11.setText(str(row_cnt + 1) + ' results')
        QtGui.QApplication.processEvents()



    def sort_table(self, column):

        fnx_mn = self.nx_min.value()
        fnx_mx = self.nx_max.value()
        fny_mn = self.ny_min.value()
        fny_mx = self.ny_max.value()
        fnz_mn = self.nz_min.value()
        fnz_mx = self.nz_max.value()
        fnpex_mn = self.npex_min.value()
        fnpex_mx = self.npex_max.value()
        fnpey_mn = self.npex_min.value()
        fnpey_mx = self.npey_max.value()
        fnpxnpy_mn = self.npxnpy_min.value()
        fnpxnpy_mx = self.npxnpy_max.value()
        fnp_mn = self.np_min.value()
        fnp_mx = self.np_max.value()
        fngpts_mn = self.ngpts_min.value()
        fngpts_mx = self.ngpts_max.value()
        nxpex_mn = self.nxpex_min.value()
        nxpex_mx = self.nxpex_max.value()
        nypey_mn = self.nypey_min.value()
        nypey_mx = self.nypey_max.value()


        if column == str("nx"):
            sorted_col = "nx"

            if self.nx_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("ny"):
            sorted_col = "ny"

            if self.ny_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.nx_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nz"):
            sorted_col = "nz"

            if self.nz_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.nx_table.setChecked(False)
            self.ny_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("npex"):
            sorted_col = "npex"

            if self.npex_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("npey"):
            sorted_col = "npey"

            if self.npey_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"


            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("npexnpey"):
            sorted_col = "npxnpy"

            if self.npexnpey_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.np_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("np"):
            sorted_col = "np"

            if self.np_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("ngpts"):
            sorted_col = "ngpts"

            if self.ngpts_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.nxpex_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nxnpex"):
            sorted_col = "nxnpex"

            if self.nxpex_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nypey_table.setChecked(False)

        if column == str("nynpey"):
            sorted_col = "nynpey"

            if self.nypey_table.isChecked() is True:
                order = " DESC"

            else:

                order = " ASC"

            self.ny_table.setChecked(False)
            self.nz_table.setChecked(False)
            self.npex_table.setChecked(False)
            self.npey_table.setChecked(False)
            self.npexnpey_table.setChecked(False)
            self.np_table.setChecked(False)
            self.nx_table.setChecked(False)
            self.ngpts_table.setChecked(False)
            self.nxpex_table.setChecked(False)

        else:
            pass

        pathx = configwr.read_config()
        pathx = pathx[19]


        conn = sqlite3.connect(pathx + "/.palm_gf_data.db")
        c = conn.cursor()
        c.execute("SELECT * FROM grid_current  WHERE nx <= " + str(fnx_mx) + " AND nx >= "  + str(fnx_mn) + " AND ny <= " + str(fny_mx) + " AND ny >= " + str(fny_mn) + " AND nz <= " + str(fnz_mx) +
                  " AND nz >= " + str(fnz_mn) + " AND npex <= " + str(fnpex_mx) + " AND npex >= " +
        str(fnpex_mn) + " AND npey <= " + str(fnpey_mx) + " AND npey >= " + str(fnpey_mn) + " AND "
                  "npxnpy <= " + str(fnpxnpy_mx) + " AND npxnpy >= " + str(fnpxnpy_mn) + " AND np <= " + str(fnp_mx) + " AND np >= " + str(fnp_mn) + " AND ngpts <= " + str(fngpts_mx) + " AND ngpts >= " + str(fngpts_mn) +
        " AND nxnpex <= " + str(nxpex_mx) + " AND nxnpex >= " + str(nxpex_mn) + " AND nynpey <= " + str(nypey_mx) + " AND nynpey >= " + str(nypey_mn) +
        " ORDER BY " + str(sorted_col) + str(order))



        sorted = c.fetchall()


        c.close()
        conn.close()
        self.tableWidget.setRowCount(len(sorted))

        for row_indx in range(0,len(sorted)):

            for col_indx in range(0,10):
                row = sorted[row_indx]
                value = row[col_indx]
                if col_indx == 7:
                    self.tableWidget.setItem(row_indx, col_indx, QtGui.QTableWidgetItem(str("%.1e" % value)))
                else:
                    self.tableWidget.setItem(row_indx, col_indx, QtGui.QTableWidgetItem(str(value)))



    def sort_order(self):

        sorted = 0
        fnx_mn = self.nx_min.value()
        fnx_mx = self.nx_max.value()
        fny_mn = self.ny_min.value()
        fny_mx = self.ny_max.value()
        fnz_mn = self.nz_min.value()
        fnz_mx = self.nz_max.value()
        fnpex_mn = self.npex_min.value()
        fnpex_mx = self.npex_max.value()
        fnpey_mn = self.npex_min.value()
        fnpey_mx = self.npey_max.value()
        fnpxnpy_mn = self.npxnpy_min.value()
        fnpxnpy_mx = self.npxnpy_max.value()
        fnp_mn = self.np_min.value()
        fnp_mx = self.np_max.value()
        fngpts_mn = self.ngpts_min.value()
        fngpts_mx = self.ngpts_max.value()

        if str(self.Sortvariable_box.currentIndex()) == str(1):


            sorted_col = "nx"
        if str(self.Sortvariable_box.currentIndex()) == str(2):

            sorted_col = "ny"
        if str(self.Sortvariable_box.currentIndex()) == str(3):

            sorted_col = "nz"
        if str(self.Sortvariable_box.currentIndex()) == str(4):

            sorted_col = "npex"
        if str(self.Sortvariable_box.currentIndex()) == str(5):

            sorted_col = "npey"
        if str(self.Sortvariable_box.currentIndex()) == str(6):

            sorted_col = "npxnpy"
        if str(self.Sortvariable_box.currentIndex()) == str(7):

            sorted_col = "np"
        if str(self.Sortvariable_box.currentIndex()) == str(8):

            sorted_col = "ngpts"

        #print(self.Sortvariable_box.currentIndex())

        if str(self.Sortorder_box.currentIndex()) == str(1):

            order = " ASC"

        if str(self.Sortorder_box.currentIndex()) == str(2):

            order = " DESC"



        conn = sqlite3.connect("/localdata/.palm_gf_data.db")
        c = conn.cursor()
        c.execute("SELECT * FROM grid_current  WHERE nx <= " + str(fnx_mx) + " AND nx >= "  + str(fnx_mn) + " AND ny <= " + str(fny_mx) + " AND ny >= " + str(fny_mn) + " AND nz <= " + str(fnz_mx) +
                  " AND nz >= " + str(fnz_mn) + " AND npex <= " + str(fnpex_mx) + " AND npex >= " +
        str(fnpex_mn) + " AND npey <= " + str(fnpey_mx) + " AND npey >= " + str(fnpey_mn) + " AND "
                  "npxnpy <= " + str(fnpxnpy_mx) + " AND npxnpy >= " + str(fnpxnpy_mn) + " AND np <= " + str(fnp_mx) + " AND np >= " + str(fnp_mn) + " AND ngpts <= " + str(fngpts_mx) + " AND ngpts >= " + str(fngpts_mn) +
        " ORDER BY " + str(sorted_col) + str(order))



        sorted = c.fetchall()

        c.close()
        conn.close()
        self.tableWidget.setRowCount(len(sorted))

        for row_indx in range(0,len(sorted)):

            for col_indx in range(0,8):
                row = sorted[row_indx]
                value = row[col_indx]
                self.tableWidget.setItem(row_indx, col_indx, QtGui.QTableWidgetItem(str(var)))






        #print(len(sorted))
        #print(sorted)

    def get_path(wildcard):
        import wx
        app = wx.App(None)
        style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        dialog = wx.FileDialog(None, 'Open', wildcard='*.db')
        if dialog.ShowModal() == wx.ID_OK:
            path = dialog.GetPath()

            print(path, file)
        else:
            path = None
        dialog.Destroy()
        return











if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp3()
    window.setWindowTitle('Gridfilter')
    window.show()
    sys.exit(app.exec_())

