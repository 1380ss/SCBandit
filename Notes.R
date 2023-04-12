# Settings
options(max.print = 50)


#progress bar
library(progress)
# Initializes the progress bar
probar <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = B,
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)      # Width of the progress bar
for(i in 1:B) {
  probar$tick()}