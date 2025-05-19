# Minimal test view for mview

# define data frame with 1000 points
N <- 1000
df <- data.frame(
  contig = rep("a1_c1", N),
  coord = 400 * runif(N),
  value = runif(N)
)

# Register a test points profile with minimal data
points_profile(
  id = "points",
  name = "Points",
  data = df
)

axis_profile()
