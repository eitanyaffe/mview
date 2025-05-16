# Minimal test view for mview

# Register a test points profile with minimal data
points_profile(
  id = "points1",
  name = "Points 1",
  data = data.frame(
    contig = c("c1", "c1", "c1"),
    coord = c(10, 20, 30),
    value = c(1, 2, 3)
  )
)

points_profile(
  id = "points2",
  name = "Points 2",
  data = data.frame(
    contig = c("c1", "c1", "c1"),
    coord = c(15, 25, 35),
    value = c(4, 5, 6)
  )
)

axis_profile()
