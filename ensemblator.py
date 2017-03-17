#!/usr/bin/env python2

from Tkinter import *
import tkFileDialog
import tkMessageBox
from ensemblator_core import *


class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.see(END)
        self.widget.configure(state="disabled")
        self.widget.update_idletasks()


class options:
    def __init__(self):
        self.input = []
        self.outout = ""
        self.pwd = ""
        self.permissive = False
        self.semipermissive = 0
        self.align = False
        self.template = ""
        self.chain = "X"
        self.model = "0"
        self.percent = 0.7
        self.dcut = 2.5
        self.auto = False
        self.maxclust = 6
        self.groupm = StrVar()
        self.groupn = StrVar()
        self.avg = False
        self.color = False



#### GUI STUFF ####

class Utility:
    def __init__(self):
        self.applocation = os.path.dirname(sys.argv[0]);
        self.last_dir_accessed = os.path.realpath(os.environ['HOME']);
        self.dir_to_save = ""

    def select_output_dir(self,label_to_update):
        dir_to_save = tkFileDialog.askdirectory(initialdir = self.last_dir_accessed, title = "Select output directory.");
        dir_to_save_label = "..." + dir_to_save[len(dir_to_save) - 22:]

        self.last_dir_accessed = dir_to_save

        label_to_update["text"] = dir_to_save_label
        return(dir_to_save);

    def move_files(self, from_to_dict, base_path):
        for from_file in from_to_dict:
            if(os.path.isfile(os.path.join(base_path,from_file))):
                shutil.move(os.path.join(base_path,from_file), os.path.join(base_path,from_to_dict[from_file]));

    def select_input_files(self, min_number, label_to_update):
        chosenfiles = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed, multiple = True, title="Select input pdb files.");
        toRet = list();

        firstname = os.path.basename(chosenfiles[0]);
        lastname = os.path.basename(chosenfiles[len(chosenfiles)-1]);
        label_to_update["text"] = firstname + " ... " + lastname
        if self.dir_to_save != "":
            self.last_dir_accessed = os.path.dirname(self.dir_to_save)
        else:
            self.last_dir_accessed = os.path.dirname(chosenfiles[0])
        toRet = chosenfiles;

        return(toRet);

    def select_input_file(self, label_to_update):
        chosenfile = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed,
                                                  multiple = False,
                                                  title="Select template file. Dissimilar chains will be removed from the ensemble.");
        toRet = "";
        if(chosenfile != ""):
            firstname = os.path.basename(chosenfile);
            label_to_update["text"] = firstname
            print self.dir_to_save
            if self.dir_to_save != "":
                self.last_dir_accessed = (self.dir_to_save)
            else:
                self.last_dir_accessed = os.path.dirname(chosenfile)
            toRet = chosenfile;
        return(toRet);

    def select_input_ensemble(self, label_to_update):
        chosenfile = tkFileDialog.askopenfilename(initialdir = self.last_dir_accessed,
                                                  multiple = False,
                                                  title="Select input ensemble (must have been prepared using the Ensemblator.");
        toRet = "";
        if(chosenfile != ""):
            firstname = os.path.basename(chosenfile);
            label_to_update["text"] = firstname
            if self.dir_to_save != "":
                self.last_dir_accessed = (self.dir_to_save)
            else:
                self.last_dir_accessed = os.path.dirname(chosenfile)
            toRet = chosenfile;
        return(toRet);


    def prepare_input_command_run(self, options, label_to_update, window_to_update):

        try:
            prepare_input(options)
        except Exception as e:
            tkMessageBox.showerror("Oops.", str(e));

        label_to_update["text"] = "Done!";
        window_to_update.update();

    def analyze_command_run(self, options, label_to_update, window_to_update):

        try:
            analyze(options)
        except Exception as e:
            tkMessageBox.showerror("Oops.", str(e));

        label_to_update["text"] = "Done!";
        window_to_update.update();


class MainRoot:
    def __init__(self):
        self.rootWindow = Tk();
        self.rootWindow.wm_title("The Ensemblator");
        self.setup_gui();
        self.utility = Utility();
        self.rootWindow.mainloop();

    def setup_gui(self):
        self.add_row("Prepare Input",
                     "Prepare ensemble for analysis.",
                     self.prepare,
                     0
                     )
        self.add_row("Analyze",
                     "Analyze prepared ensemble.",
                     self.analyze,
                     1
                     )
        self.add_row("Exit",
                     "Exit application.",
                     self.exit,
                     2
                     )

        t1 = Text(self.rootWindow)
        sys.stdout = TextRedirector(t1, "stdout")
        sys.stderr = TextRedirector(t1, "stderr")
        t1.grid(row = 0, rowspan = 3, column = 2, sticky=E)


    def add_row(self, button_name, label_text, func_to_call, rownum):

        frame1 = Frame(self.rootWindow);
        #frame1.grid(row=0, column = 0);
        fr1_button = Button(self.rootWindow,
                            text = button_name,
                            command = func_to_call
                            )
        fr1_button.grid(row = rownum, column = 0, sticky = E+W);
        fr1_label = Label(self.rootWindow, text = label_text);
        fr1_label.grid(row = rownum, column = 1, sticky = W);

    def prepare(self):
        subwindow = Prepare(self.utility);
    def analyze(self):
        subwindow = Analyze(self.utility);
    def exit(self):
        sys.exit("Have a nice day!")



class Prepare:
    def __init__(self, utility):
        self.utility = utility;
        self.rootWindow = Tk();
        self.rootWindow.wm_title("Prepare Input");
        self.output = StringVar(self.rootWindow);
        self.output.set("prepared_ensemble.pdb");
        self.break_num = IntVar(self.rootWindow);
        self.break_num.set(0);
        self.gap_setting = StringVar(self.rootWindow);
        self.gap_setting.set("None");
        self.align = IntVar(self.rootWindow)
        self.align.set(0)
        self.template = ""
        self.chain = StringVar(self.rootWindow)
        self.chain.set("X")
        self.model = StringVar(self.rootWindow)
        self.model.set("0")
        self.percent = StringVar(self.rootWindow)
        self.percent.set("70")


        self.inputfiles = list();
        self.dir_to_save = "";
        self.setup_gui();
        self.rootWindow.mainloop();

    def setup_gui(self):

        select_button = Button(self.rootWindow,
                               text = "Select Input Files",
                               command = self.select_files
                               )
        select_button.grid(row=0,
                           column = 0,
                           columnspan = 2,
                           sticky = E+W
                           );

        files_selected_label = Label(self.rootWindow,
                                     text = "Files Selected: "
                                     )
        files_selected_label.grid(row = 1,
                                  column = 0,
                                  sticky = W
                                  )
        self.files_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.files_selected.grid(row = 1,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )


        pwd_button = Button(self.rootWindow,
                            text = "Select Working Directory",
                            command = self.select_pwd
                            )
        pwd_button.grid(row=2,
                        column = 0,
                        columnspan = 2,
                        sticky = E+W
                        )

        pwd_selected_label = Label(self.rootWindow,
                                   text = "Working Directory: "
                                   )
        pwd_selected_label.grid(row = 3,
                                column = 0,
                                sticky = W
                                )
        self.pwd_selected = Label(self.rootWindow,
                                  text = "None           "
                                  )
        self.pwd_selected.grid(row = 3,
                               column = 1,
                               columnspan = 3,
                               sticky = W
                               )


        output_label = Label(self.rootWindow,
                             text = "Ensemble output filename: "
                             )
        output_label.grid(row = 4,
                          column = 0,
                          sticky = W
                          )
        output = Entry(self.rootWindow,
                       textvariable = self.output
                       )
        output.grid(row = 4,
                    column = 1,
                    sticky=W
                    )


        gaps_label = Label(self.rootWindow,
                           text = "Chain-breaks permitted? "
                           )
        gaps_label.grid(row = 5,
                        column = 0,
                        sticky = W
                        )
        gaps_dropdown = OptionMenu(self.rootWindow,
                                   self.gap_setting,
                                   "None",
                                   "Some",
                                   "All",
                                   command=self.Semipermissive
                                   )
        gaps_dropdown.grid(row = 5,
                           column = 1,
                           sticky = E
                           )

        break_num_label = Label(self.rootWindow,
                                text = "# of chain-breaks permitted: "
                                )
        break_num_label.grid(row = 6,
                             column = 0,
                             sticky = W
                             )
        self.break_num_entry = Entry(self.rootWindow,
                                     textvariable = self.break_num
                                     )
        self.break_num_entry.grid(row = 6,
                                  column = 1,
                                  sticky=W
                                  )
        self.break_num_entry.configure(state='disable')


        align_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Perform sequence alignment",
                                        variable = self.align,
                                        command = self.align_check
                                        )
        align_checkbutton.grid(row = 7,
                               column = 0,
                               sticky = W+E
                               )




        self.template_entry = Button(self.rootWindow,
                                     text = "Select Template File",
                                     command = self.select_file
                                     )
        self.template_entry.grid(row = 8,
                                  column = 0,
                                  sticky=E
                                  )
        self.template_entry.configure(state='disable')

        self.template_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.template_selected.grid(row = 8,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )


        chain_label = Label(self.rootWindow,
                                text = "Chain ID for template: "
                                )
        chain_label.grid(row = 9,
                             column = 0,
                             sticky = W
                             )
        self.chain_entry = Entry(self.rootWindow,
                                     textvariable = self.chain
                                     )
        self.chain_entry.grid(row = 9,
                                  column = 1,
                                  sticky=W
                                  )
        self.chain_entry.configure(state='disable')




        model_label = Label(self.rootWindow,
                                text = "Model ID for template: "
                                )
        model_label.grid(row = 10,
                             column = 0,
                             sticky = W
                             )
        self.model_entry = Entry(self.rootWindow,
                                 text = "0",
                                 textvariable = self.model
                                     )
        self.model_entry.grid(row = 10,
                                  column = 1,
                                  sticky=W
                                  )
        self.model_entry.configure(state='disable')

        percent_label = Label(self.rootWindow,
                                text = "Percent ID Cutoff: "
                                )
        percent_label.grid(row = 11,
                             column = 0,
                             sticky = W
                             )
        self.percent_entry = Entry(self.rootWindow,
                                   text = "0.7",
                                   textvariable = self.percent
                                   )
        self.percent_entry.grid(row = 11,
                                  column = 1,
                                  sticky=W
                                  )
        self.percent_entry.configure(state='disable')



        go_button = Button(self.rootWindow, text = "Go!", command = self.execute);
        go_button.grid(row=12, column = 0, columnspan = 2, sticky = E+W);

        status_label = Label(self.rootWindow, text = "Status: ");
        status_label.grid(row = 13, column = 0, sticky = W);
        self.status = Label(self.rootWindow, text = "Waiting");
        self.status.grid(row = 13, column = 1, columnspan = 3, sticky = W);


    def align_check(self):
        if self.align.get() == 1:
            self.template_entry.configure(state = 'normal')
            self.chain_entry.configure(state = 'normal')
            self.model_entry.configure(state = 'normal')
            self.percent_entry.configure(state = 'normal')
        else:
            self.template_entry.configure(state = 'disable')
            self.chain_entry.configure(state = 'disable')
            self.model_entry.configure(state = 'disable')
            self.percent_entry.configure(state = 'disable')

    def Semipermissive(self, event):
        if self.gap_setting.get() == "Some":
            self.break_num_entry.configure(state='normal')
        else:
            self.break_num_entry.configure(state='disable')


    def select_files(self):
        self.status["text"] = "Waiting";
        self.inputfiles = self.utility.select_input_files(1,self.files_selected)

    def select_file(self):
        self.status["text"] = "Waiting";
        self.template = self.utility.select_input_file(self.template_selected)

    def select_pwd(self):
        self.status["text"] = "Waiting";
        self.dir_to_save = self.utility.select_output_dir(self.pwd_selected)

    def execute(self):
        self.status["text"] = "Waiting";

        break_num = int(self.break_num.get());
        gap_setting = self.gap_setting.get();



        ## Manually specify the arguments!

        # shared arguments
        options.pwd = self.dir_to_save

        # prepare input arguments
        options.input = self.inputfiles
        options.output = str(self.output.get())
        if gap_setting == "None":
            options.permissive = False
            options.semipermissive = 0
        elif gap_setting == "All":
            options.permissive = True
            options.semipermissive = 0
        elif gap_setting == "Some":
            options.permissive = False
            options.semipermissive = break_num
        else:
            options.permissive = False
            options.semipermissive = 0


        if self.align.get() == 1:
            options.align = True
            options.template = os.path.basename(self.template)\
                [0:(len(os.path.basename(self.template)) - 4)]

            options.chain = self.chain.get()
            options.model = self.model.get()
            options.percent = float(self.percent.get())

        else:
            options.align = False



        if(self.dir_to_save != "" and len(self.inputfiles) > 0):
            if (self.align.get() == 1 and self.template != "") \
                or (self.align.get() == 0):
                self.utility.prepare_input_command_run(options,
                                                       self.status,
                                                       self.rootWindow
                                                       )
                self.rootWindow.destroy()
            else:
                print("\n\nIf you have chosen to do an alignment, you must"
                      " also select a template file to use as your guide."
                      " Models that are less similar than the percent identity"
                      " cutoff will not be included in the ensemble."
                      )


        else:
            print("\n\nPlease select a working directory (to save results in)"
                  ", and some input files to make an ensemble out of!"
                  )


class Analyze:

    def __init__(self, utility):
        self.utility = utility
        self.rootWindow = Tk()
        self.rootWindow.wm_title("Analyze Ensemble")
        self.dcut = StringVar(self.rootWindow)
        self.dcut.set("2.5")
        self.dcutAuto = IntVar(self.rootWindow)
        self.dcutAuto.set(0)
        self.ensemble = ""
        self.groupm = StringVar(self.rootWindow)
        self.groupm.set("")
        self.groupn = StringVar(self.rootWindow)
        self.groupn.set("")
        self.maxclust = IntVar(self.rootWindow)
        self.maxclust.set(3)
        self.auto = IntVar(self.rootWindow)
        self.auto.set(0)
        self.color = IntVar(self.rootWindow)
        self.color.set(0)
        self.dir_to_save = ""
        self.setup_gui()
        self.rootWindow.mainloop()

    def setup_gui(self):

        pwd_button = Button(self.rootWindow,
                            text = "Select Working Directory",
                            command = self.select_pwd
                            )
        pwd_button.grid(row=0,
                        column = 0,
                        columnspan = 2,
                        sticky = E+W
                        )

        pwd_selected_label = Label(self.rootWindow,
                                   text = "Working Directory: "
                                   )
        pwd_selected_label.grid(row = 1,
                                column = 0,
                                sticky = W
                                )
        self.pwd_selected = Label(self.rootWindow,
                                  text = "None           "
                                  )
        self.pwd_selected.grid(row = 1,
                               column = 1,
                               columnspan = 3,
                               sticky = W
                               )



        select_button = Button(self.rootWindow,
                               text = "Select Input Ensemble",
                               command = self.select_ensemble
                               )
        select_button.grid(row=2,
                           column = 0,
                           columnspan = 2,
                           sticky = E+W
                           );

        file_selected_label = Label(self.rootWindow,
                                     text = "Ensemble Selected: "
                                     )
        file_selected_label.grid(row = 3,
                                  column = 0,
                                  sticky = W
                                  )
        self.ensemble_selected = Label(self.rootWindow,
                                    text = "None           "
                                    )
        self.ensemble_selected.grid(row = 3,
                                 column = 1,
                                 columnspan = 3,
                                 sticky = W
                                 )



        dcut_label = Label(self.rootWindow,
                             text = "Core cutoff distance: "
                             )
        dcut_label.grid(row = 5,
                          column = 0,
                          sticky = W
                          )
        self.dcut_entry = Entry(self.rootWindow,
                       textvariable = self.dcut
                       )
        self.dcut_entry.grid(row = 5,
                    column = 1,
                    sticky=W
                    )

        dcut_auto_checkbutton = Checkbutton(self.rootWindow,
                                       text="Automatic Cutoff Distance Search",
                                       variable=self.dcutAuto,
                                       command=self.dcut_auto_check
                                       )
        dcut_auto_checkbutton.grid(row=4,
                              column=0,
                              columnspan=2,
                              sticky=E
                              )



        groupm_label = Label(self.rootWindow,
                                text = "Group M models: "
                                )
        groupm_label.grid(row = 6,
                             column = 0,
                             sticky = W
                             )
        self.groupm_entry = Entry(self.rootWindow,
                                     textvariable = self.groupm
                                     )
        self.groupm_entry.grid(row = 6,
                                  column = 1,
                                  sticky=W
                                  )


        groupn_label = Label(self.rootWindow,
                                text = "Group N models: "
                                )
        groupn_label.grid(row = 7,
                             column = 0,
                             sticky = W
                             )
        self.groupn_entry = Entry(self.rootWindow,
                                 text = "0",
                                 textvariable = self.groupn
                                     )
        self.groupn_entry.grid(row = 7,
                                  column = 1,
                                  sticky=W
                                  )






        auto_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Perform automatic clustering",
                                        variable = self.auto,
                                        command = self.auto_check
                                        )
        auto_checkbutton.grid(row = 8,
                               column = 0,
                               columnspan = 2,
                               sticky = W
                               )



        maxclust_label = Label(self.rootWindow,
                                text = "Max # of clusters to search for: "
                                )
        maxclust_label.grid(row = 9,
                             column = 0,
                             sticky = W
                             )
        self.maxclust_entry = Entry(self.rootWindow,
                                     textvariable = self.maxclust
                                     )
        self.maxclust_entry.grid(row = 9,
                                  column = 1,
                                  sticky=W
                                  )
        self.maxclust_entry.configure(state='disable')


        color_checkbutton = Checkbutton(self.rootWindow,
                                        text = "Set b-factors in final " +\
                                               "ensemble equal to inter-LODR" +\
                                               " (or group M LODR).",
                                        variable = self.color
                                        )
        color_checkbutton.grid(row = 11,
                               column = 0,
                               columnspan = 2,
                               sticky = W
                               )




        go_button_analyze = Button(self.rootWindow,
                                   text = "Analyze!",
                                   command = self.execute_analyze
                                   )
        go_button_analyze.grid(row=12,
                               column = 0,
                               columnspan = 2,
                               sticky = E+W
                               )

        status_label = Label(self.rootWindow, text = "Status: ");
        status_label.grid(row = 13, column = 0, sticky = W);
        self.status = Label(self.rootWindow, text = "Waiting");
        self.status.grid(row = 13, column = 1, columnspan = 3, sticky = W);

    def auto_check(self):
        if self.auto.get() == 1:
            self.groupm_entry.configure(state = 'disable')
            self.groupn_entry.configure(state = 'disable')
            self.maxclust_entry.configure(state='normal')
        else:
            self.groupm_entry.configure(state = 'normal')
            self.groupn_entry.configure(state = 'normal')
            self.maxclust_entry.configure(state='disable')


    def dcut_auto_check(self):
        if self.dcutAuto.get() == 1:
            self.dcut_entry.configure(state = 'disable')
        else:
            self.dcut_entry.configure(state = 'normal')


    def select_ensemble(self):
        self.status["text"] = "Waiting";
        self.ensemble = self.utility.select_input_ensemble(self.ensemble_selected)

    def select_pwd(self):
        self.status["text"] = "Waiting";
        self.dir_to_save = self.utility.select_output_dir(self.pwd_selected)

    # parser for group numbers
    def parse_range(self, value):
        result = set()
        for part in value.split(','):
            x = part.split('-')
            result.update(range(int(x[0]), int(x[-1]) + 1))
        return(sorted(result))


    def execute_analyze(self):
        self.status["text"] = "Waiting";

        ## Manually specify the arguments!
        options.pwd = self.dir_to_save


        # analyze arguments
        options.input = self.ensemble
        options.dcut = float(self.dcut.get())
        options.dcutAuto = self.dcutAuto.get()

        options.auto = self.auto.get()
        if options.auto == 1:
            options.maxclust = self.maxclust.get()
        else:
            try:
                options.groupm = self.parse_range(str(self.groupm.get()))
            except:
                print "\n\nPlease select a group M, or choose the auto option."
                print "\nGroups should be formatted like so: 0,3,5,6-9,13\n\n"

            try:
                options.groupn = self.parse_range(str(self.groupn.get()))
            except:
                pass

        options.color = True if (self.color.get() == 1) else False



        if(self.dir_to_save != "" and self.ensemble != "") and\
           (self.auto.get() == 1 or self.groupm.get() != ""):


                self.utility.analyze_command_run(options,
                                                       self.status,
                                                       self.rootWindow
                                                       )
                self.rootWindow.destroy()
        else:
            print("\n\nPlease select a working directory (to save results in)"
                  ", a prepared ensemble to analyze, and either select the "
                  "auto-cluster option, or define at least a group M to "
                  "analyze. Then click the button again!"
                  )




main_window = MainRoot();
