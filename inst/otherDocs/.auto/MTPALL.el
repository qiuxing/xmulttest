(TeX-add-style-hook "MTPALL"
 (lambda ()
    (LaTeX-add-bibliographies
     "xmulttest")
    (LaTeX-add-environments
     "theorem"
     "procedure")
    (LaTeX-add-labels
     "f:cytoPlot"
     "f:cytoMargPlot"
     "f:cytogfwer"
     "f:cytotppfp"
     "f:cytofdr"
     "f:mbPlot"
     "f:coxphPlot")
    (TeX-add-symbols
     '("Rclass" 1)
     '("Robject" 1)
     '("Rpackage" 1)
     "RR"
     "ZZ"
     "NN")
    (TeX-run-style-hooks
     "natbib"
     "authoryear"
     "round"
     "comment"
     "color"
     "amsmath"
     "hyperref"
     "amsfonts"
     "Sweave"
     "graphicx"
     "latex2e"
     "art11"
     "article"
     "11pt")))

