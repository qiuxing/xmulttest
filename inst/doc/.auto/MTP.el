(TeX-add-style-hook "MTP"
 (lambda ()
    (LaTeX-add-bibliographies
     "xmulttest")
    (LaTeX-add-environments
     "theorem"
     "procedure")
    (LaTeX-add-labels
     "anal:mult:xmulttest"
     "anal:mult:s:intro"
     "anal:mult:s:methods"
     "anal:mult:s:framework"
     "anal:mult:e:tstat"
     "anal:mult:t:TypeIandII"
     "anal:mult:e:gFWER"
     "anal:mult:e:PCER"
     "anal:mult:e:TPPFP"
     "anal:mult:e:FDR"
     "anal:mult:s:nullDistn"
     "anal:mult:proc:boot"
     "anal:mult:s:SS"
     "anal:mult:e:SScut"
     "anal:mult:e:SSquant"
     "anal:mult:s:SD"
     "anal:mult:e:SDmaxT"
     "anal:mult:e:SDminP"
     "anal:mult:s:AMTP"
     "anal:mult:e:adjpgFWER"
     "anal:mult:e:augTPPFP"
     "anal:mult:e:adjpTPPFP"
     "anal:mult:s:software"
     "anal:mult:s:MTP"
     "anal:mult:s:summaries"
     "anal:mult:s:design"
     "anal:mult:s:disc")
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

