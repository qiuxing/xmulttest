(TeX-add-style-hook "multtest"
 (lambda ()
    (LaTeX-add-bibliographies
     "xmulttest")
    (LaTeX-add-labels
     "fig:mtQQ"
     "fig:mtNumDen"
     "fig:mtpvsr"
     "fig:mtpvst")
    (TeX-run-style-hooks
     "hyperref"
     "natbib"
     "authoryear"
     "round"
     "graphicx"
     "fullpage"
     "epsfig"
     "amsmath"
     "latex2e"
     "art11"
     "article"
     "11pt")))

