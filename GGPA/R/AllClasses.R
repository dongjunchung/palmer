
# GGPA class definition

setClass( Class="GGPA",
  representation=representation(
    fit="list",
    summary="list",
	  setting="list",
    gwasPval="matrix"
  )
)
