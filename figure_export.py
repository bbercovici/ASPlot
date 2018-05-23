import os


# I need to generate the correct calls to 
# \includegraphics from here




'''
Inputs:
-------
- inputfolder : path to folder enclosing all saved 
'''
def generate_latex_skeleton(inputfolder):


	latex_calls = []

    for folder in os.walk(inputfolder) :
        for subfolder in folder[1]:
            foldername = folder[0] + subfolder + "/"

            latex_calls += "\subsection{ Case " + subfolder + "}\n"
            latex_calls += ["\includegraphics[width = \textwidth]{" + foldername + "positions.pdf" + "}\n"]
            latex_calls += ["\includegraphics[width = \textwidth]{" + foldername + "velocities.pdf" + "}\n"]


    text_file = open("test.tex", "w")
    for line in latex_calls:
		text_file.write(line)
	text_file.close()





