.onAttach <- function(...) {
  packageStartupMessage("drord: Doubly robust estimators for ordinal outcomes")
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("drord")$Version
  )
}
