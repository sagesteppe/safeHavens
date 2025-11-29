## Append-mode globals.R updater

# (Assumes your working directory is package root, and that R/globals.R exists or not.)

if (! requireNamespace("checkhelper", quietly = TRUE)) {
  stop("Please install checkhelper first: install.packages('checkhelper')")
}
library(checkhelper)

# 1. Run checkhelper to detect missing globals
globals <- get_no_visible(path = ".", quiet = TRUE)
new_globals <- checkhelper::print_globals(globals, what = "globalVariables", show = FALSE)
# new_globals is a character vector: the lines of the globalVariables(...) call

# 2. Try to read existing R/globals.R (if present), parse existing globals
existing <- character()
globals_file <- file.path("..", "R", "globals.R")
if (file.exists(globals_file)) {
  content <- readLines(globals_file, warn = FALSE)
  # crude parse: look for utils::globalVariables or globalVariables call
  gv_lines <- grep("globalVariables\\(", content, value = TRUE)
  if (length(gv_lines) > 0) {
    # extract all quoted names inside the call(s)
    # This regex: match "ABC", possibly multiple, inside globalVariables(...).
    all_quotes <- unlist(regmatches(gv_lines,
                                   gregexpr('"(?:\\\\.|[^"\\\\])*"', gv_lines)))
    existing <- unique(gsub('"', "", all_quotes))
  }
}

# 3. Extract new globals from new_globals output
# new_globals likely contains a call like:
# globalVariables(unique(c("foo", "bar", ...)))
# So we can extract all quoted names similarly
new_quotes <- unlist(regmatches(new_globals,
                               gregexpr('"(?:\\\\.|[^"\\\\])*"', new_globals)))
new_vars <- unique(gsub('"', "", new_quotes))

# 4. Merge
merged <- sort(unique(c(existing, new_vars)))

# 5. Create the merged globalVariables call
out <- c(
  "# This file was auto-generated (or updated) by checkhelper",
  "if(getRversion() >= \"2.15.1\") {",
  "  utils::globalVariables(",
  sprintf("    c(%s)", paste0(sprintf('"%s"', merged), collapse = ", ")),
  "  )",
  "}"
)

# 6. Write to R/globals.R
cat(paste0(out, collapse = "\n"), file = globals_file)
message("Wrote ", length(merged), " global variable(s) to ", globals_file, ".")

