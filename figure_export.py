# MIT License

# Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import os


# The cleanest way to get this to work is probably to type something along the lines of 
# rsync -r -R bebe0705@fortuna.colorado.edu:/home/anfr8485/FALCON/Filter/test/advspc_results_1/*/*.pdf .
# and call generate_latex_skeleton with inputfolder == "./home/anfr8485/FALCON/Filter/test/advspc_results_1/"
# yes, you read this path right (rsync preserves the complete directory structure)

# example: generate_latex_skeleton("/Users/bbercovici/GDrive/CUBoulder/Research/code/AS/home/anfr8485/FALCON/Filter/test/advspc_results_1","test.tex")


'''
Inputs:
-------
- inputfolder : path to folder enclosing all saved cases
- outputpath : path to output tex file (ex : "figures.tex") 
'''
def generate_latex_skeleton(inputfolder,outputpath):


    latex_calls = [r"\section{Results}"]
    latex_calls += "\n"

    for folder in os.walk(inputfolder) :
        for subfolder in folder[1]:
            foldername = folder[0] + subfolder + "/"

            subsectionName = subfolder.replace('_', '\\_')

            latex_calls += r"\subsection{$\mathrm{" + subsectionName + r"}$}"
            latex_calls += "\n"
            latex_calls += "\\centering"
            latex_calls += "\n"
            latex_calls += [r"\includegraphics[width = 0.7\textwidth]{" + foldername + r"positions.pdf" + r"}"]
            latex_calls += "\n"
            latex_calls += "\\newline"
            latex_calls += "\n"
            latex_calls += "\\centering"
            latex_calls += "\n"
            latex_calls += [r"\includegraphics[width = 0.7\textwidth]{" + foldername + r"velocities.pdf" + r"}"]
            latex_calls += "\n"

            latex_calls += "\\begin{table}[!htb]"
            latex_calls += "\n"
            latex_calls += "\\begin{center}"
            latex_calls += "\n"
            latex_calls += r"\captionsetup{width=.8\linewidth}"
            latex_calls += "\n"
            latex_calls += r"\caption{State sigmas mapped to maneuver times in km and km/s}"
            latex_calls += "\n"
            latex_calls += r"\small"
            latex_calls += "\n"
            latex_calls += "\\begin{tabular}{l c c c c c c}"
            latex_calls += "\n"
            latex_calls += "\\toprule"
            latex_calls += "\n"
            latex_calls += r"Maneuver & $\sigma_{R}$ & $\sigma_{T}$ & $\sigma_{N}$ & $\sigma_{\dot{R}}$ & $\sigma_{\dot{T}}$ & $\sigma_{\dot{N}}$ \\"
            latex_calls += "\n"
            latex_calls += r"\otoprule"
            latex_calls += "\n" 

            with open(foldername + r"sigmas.map") as f:
                lines = f.readlines()

            for line in lines:
                latex_calls += line
                latex_calls += "\n"


            latex_calls += "\\bottomrule"
            latex_calls += "\n"
            latex_calls += r"\end{tabular}"
            latex_calls += "\n"
            latex_calls += r"\end{center}"
            latex_calls += "\n"
            latex_calls += r"\end{table}"
            latex_calls += "\n"


            latex_calls += r"\clearpage"
            latex_calls += "\n"

    text_file = open(outputpath, "w")
    for line in latex_calls:
        text_file.write(line)
    text_file.close()


workingDir = 'home/anfr8485/FALCON/Filter/test/advspc_results/'
pdfDir = '/Users/afrench/Desktop/Adv Spc Report/results.tex'
generate_latex_skeleton(os.path.join(os.getcwd(), workingDir), pdfDir)






