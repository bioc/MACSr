#' @import basilisk
env_macs <- BasiliskEnvironment("env_macs", pkgname="MACSr",
                                packages = c("python=3.12"),
                                pip = c("macs3==3.1.0"))
