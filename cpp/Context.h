#ifndef CONTEXT_H
#define CONTEXT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <numeric>

// Input structures
struct SegmentRow {
  std::string segment_id;
  std::string contig;
  int64_t start;
  int64_t end;
};

struct ContigRow {
  std::string contig_id;
  int64_t length;
};

// Output structures for point transformations
struct PointRow {
  std::string contig;
  int64_t coord;
  double vcoord;  // for contig2view results
  int input_index;  // 1-based index of input row that produced this output
};

// Output structures for interval transformations
struct IntervalRow {
  std::string contig;
  int64_t start;
  int64_t end;
  double vstart;  // for contig2view results
  double vend;
  bool trim_left;
  bool trim_right;
  int n_segments;
  std::string segment_ids;  // comma-separated when merged
  int input_index;  // 1-based index of input row that produced this output
};

// PlottedSegment structure (used in public interface)
struct PlottedSegment {
  std::string segment_id;
  std::string contig;
  int64_t contig_start;  // original start on contig
  int64_t contig_end;    // original end on contig
  int64_t vcoord_start;  // cumulative start in vcoord space
  int64_t vcoord_end;    // cumulative end in vcoord space
};

class Context {
public:
  // Default constructor (does nothing)
  Context();
  
  // Initialize/reinitialize tables
  void set_tables(const std::vector<SegmentRow>& segment_table, 
                  const std::vector<ContigRow>& contig_table);
  
  // Update with user-selected segments
  void update_selected_segments(const std::vector<SegmentRow>& selected_segments, 
                                const std::vector<double>& zoom = std::vector<double>());
  
  // Point transformations
  std::vector<PointRow> contig2view_point(const std::vector<PointRow>& input);
  std::vector<PointRow> view2contig_point(const std::vector<PointRow>& input);
  
  // Interval transformations
  std::vector<IntervalRow> contig2view_interval(const std::vector<IntervalRow>& input, 
                                                 const std::string& crossing_intervals = "drop");
  std::vector<IntervalRow> view2contig_interval(const std::vector<IntervalRow>& input);
  
  // Get view intervals with optional xlim clipping and merging
  std::vector<IntervalRow> get_view_intervals(bool limit_to_zoom = true, bool merge_adjacent = true);
  
  // Filtering methods
  std::vector<PointRow> filter_coords(const std::vector<PointRow>& input, 
                                      const std::vector<double>& xlim = std::vector<double>());
  std::vector<IntervalRow> filter_intervals(const std::vector<IntervalRow>& input, 
                                             const std::vector<double>& xlim = std::vector<double>(),
                                             bool merge_adjacent = false);
  
  // Check if coords are in view (returns boolean for each input)
  std::vector<bool> coords_in_view(const std::vector<PointRow>& input, bool limit_to_zoom = false);
  
  // Accessor methods
  std::vector<PlottedSegment> get_plotted_segments();
  std::vector<double> get_xlim();
  std::vector<std::string> get_segment_ids(bool limit_to_zoom = false);
  
private:
  // Internal data structures
  struct Segment {
    std::string segment_id;
    std::string contig;
    int64_t start;  // 1-based start on contig
    int64_t end;    // 1-based end on contig (inclusive)
    int64_t length;
    int64_t acoord_start;  // cumulative start in acoord space
  };
  
  struct Contig {
    std::string contig_id;
    int64_t length;
    int64_t acoord_start;  // cumulative start in acoord space
  };
  
  // Fixed reference tables
  // NOTE: Coordinate system design
  // - acoord (assembly coordinate): 0-based global coordinate system shared by both contigs and segments
  // - Both contigs and segments map directly to acoord independently (not hierarchical for conversions)
  // - While segments belong to contigs structurally, for coordinate conversions we treat them as
  //   two separate ways of partitioning the same acoord space
  // - vcoord (view coordinate): 0-based coordinate spanning only plotted/selected segments
  // - To convert contig <-> vcoord, we must go through segments because vcoord only includes plotted segments
  std::vector<Segment> full_segments_;
  std::vector<Contig> contigs_;
  std::unordered_map<std::string, size_t> segment_id_to_index_;
  std::unordered_map<std::string, size_t> contig_id_to_index_;
  
  // Prefix sums for binary search (exclusive ends: cs[i] = cumulative length up to element i)
  std::vector<int64_t> all_contig_cs;
  std::vector<int64_t> all_segment_cs;
  
  // Current view state
  std::vector<PlottedSegment> plotted_segments_;
  std::vector<double> xlim_;
  std::unordered_map<std::string, size_t> plotted_segment_id_to_index_;
  std::vector<int64_t> plotted_segment_cs;
  
  // Helper methods
  int64_t contig2acoord(const std::string& contig, int64_t coord);
  std::pair<std::string, int64_t> acoord2contig(int64_t acoord);
  std::pair<std::string, int64_t> acoord2segment(int64_t acoord);
  int64_t segment2acoord(const std::string& segment, int64_t coord);
  int64_t segment2vcoord(const std::string& segment, int64_t coord);
  std::pair<std::string, int64_t> vcoord2segment(int64_t vcoord);
  
  void build_full_segments(const std::vector<SegmentRow>& segment_table);
  void build_contigs(const std::vector<ContigRow>& contig_table);
  void update_plotted_segments(const std::vector<SegmentRow>& selected_segments);
  void update_xlim(const std::vector<double>& zoom);
};

#endif // CONTEXT_H

