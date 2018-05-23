import os

'''
Inputs:
-------
- inputfolder : path to folder enclosing all saved cases
- outputpath : path to output tex file (ex : "figures.tex") 
'''
def generate_latex_skeleton(inputfolder,outputpath):


    latex_calls = []

    for folder in os.walk(inputfolder) :
        for subfolder in folder[1]:
            foldername = folder[0] + subfolder + "/"

            latex_calls += r"\subsection{$\mathrm{" + subfolder + r"}}\n"
            latex_calls += [r"\includegraphics[width = \textwidth]{" + foldername + r"positions.pdf" + r"}\n"]
            latex_calls += [r"\includegraphics[width = \textwidth]{" + foldername + r"velocities.pdf" + r"}\n"]


    text_file = open(outputpath, "w")
    for line in latex_calls:
        text_file.write(line)
    text_file.close()





