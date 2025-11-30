// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include "Context.h"
#include <stdexcept>

using namespace Rcpp;

// Context class wrapper using XPtr
typedef Rcpp::XPtr<Context> ContextPtr;

namespace {
Context* require_context(ContextPtr ctx, const char* func_name) {
  if (ctx.isNULL() || ctx.get() == nullptr) {
    Rcpp::stop("Context pointer is NULL in %s. Call context_create() before using this function.", func_name);
  }
  return ctx.get();
}

std::vector<double> to_std_vector(NumericVector vec) {
  std::vector<double> out;
  for (int i = 0; i < vec.size(); ++i) {
    out.push_back(vec[i]);
  }
  return out;
}

// DataFrame to STL conversion functions
std::vector<SegmentRow> dataframe_to_segment_rows(DataFrame df) {
  CharacterVector segment_ids = df["segment"];
  CharacterVector contigs = df["contig"];
  IntegerVector starts = df["start"];
  IntegerVector ends = df["end"];
  
  int n = df.nrows();
  std::vector<SegmentRow> result;
  result.reserve(n);
  
  for (int i = 0; i < n; i++) {
    SegmentRow row;
    row.segment_id = as<std::string>(segment_ids[i]);
    row.contig = as<std::string>(contigs[i]);
    row.start = starts[i];
    row.end = ends[i];
    result.push_back(row);
  }
  
  return result;
}

std::vector<ContigRow> dataframe_to_contig_rows(DataFrame df) {
  CharacterVector contig_ids = df["contig"];
  IntegerVector lengths = df["length"];
  
  int n = df.nrows();
  std::vector<ContigRow> result;
  result.reserve(n);
  
  for (int i = 0; i < n; i++) {
    ContigRow row;
    row.contig_id = as<std::string>(contig_ids[i]);
    row.length = lengths[i];
    result.push_back(row);
  }
  
  return result;
}

std::vector<PointRow> dataframe_to_point_rows(DataFrame df, bool has_vcoord = false) {
  int n = df.nrows();
  std::vector<PointRow> result;
  result.reserve(n);
  
  if (has_vcoord) {
    // For view2contig: only need vcoord, contig and coord will be filled by C++
    if (!df.containsElementNamed("vcoord")) {
      Rcpp::stop("Error in dataframe_to_point_rows: missing required column 'vcoord'");
    }
    NumericVector vcoords = df["vcoord"];
    for (int i = 0; i < n; i++) {
      PointRow row;
      row.contig = "";  // will be filled by view2contig
      row.coord = 0;  // will be filled by view2contig
      row.vcoord = vcoords[i];
      result.push_back(row);
    }
  } else {
    // For contig2view: need contig and coord
    if (!df.containsElementNamed("contig")) {
      Rcpp::stop("Error in dataframe_to_point_rows: missing required column 'contig'");
    }
    if (!df.containsElementNamed("coord")) {
      Rcpp::stop("Error in dataframe_to_point_rows: missing required column 'coord'");
    }
    CharacterVector contigs = df["contig"];
    IntegerVector coords = df["coord"];
    for (int i = 0; i < n; i++) {
      PointRow row;
      row.contig = as<std::string>(contigs[i]);
      row.coord = coords[i];  // R 1-based coordinates
      row.vcoord = 0.0;  // will be filled by contig2view
      result.push_back(row);
    }
  }
  
  return result;
}

std::vector<IntervalRow> dataframe_to_interval_rows(DataFrame df, bool has_vcoords = false) {
  int n = df.nrows();
  std::vector<IntervalRow> result;
  result.reserve(n);
  
  if (has_vcoords) {
    // For view2contig: only need vstart and vend, contig will be filled by C++
    if (!df.containsElementNamed("vstart")) {
      Rcpp::stop("Error in dataframe_to_interval_rows: missing required column 'vstart'");
    }
    if (!df.containsElementNamed("vend")) {
      Rcpp::stop("Error in dataframe_to_interval_rows: missing required column 'vend'");
    }
    NumericVector vstarts = df["vstart"];
    NumericVector vends = df["vend"];
    for (int i = 0; i < n; i++) {
      IntervalRow row;
      row.contig = "";  // will be filled by view2contig
      row.start = 0;  // will be filled by view2contig
      row.end = 0;
      row.vstart = vstarts[i];
      row.vend = vends[i];
      row.trim_left = false;
      row.trim_right = false;
      row.n_segments = 0;
      result.push_back(row);
    }
  } else {
    // For contig2view: need contig, start, end
    if (!df.containsElementNamed("contig")) {
      Rcpp::stop("Error in dataframe_to_interval_rows: missing required column 'contig'");
    }
    if (!df.containsElementNamed("start")) {
      Rcpp::stop("Error in dataframe_to_interval_rows: missing required column 'start'");
    }
    if (!df.containsElementNamed("end")) {
      Rcpp::stop("Error in dataframe_to_interval_rows: missing required column 'end'");
    }
    CharacterVector contigs = df["contig"];
    IntegerVector starts = df["start"];
    IntegerVector ends = df["end"];
    for (int i = 0; i < n; i++) {
      IntervalRow row;
      row.contig = as<std::string>(contigs[i]);
      row.start = starts[i];  // R 1-based coordinates
      row.end = ends[i];  // R 1-based coordinates
      row.vstart = 0.0;  // will be filled by contig2view
      row.vend = 0.0;
      row.trim_left = false;
      row.trim_right = false;
      row.n_segments = 0;
      result.push_back(row);
    }
  }
  
  return result;
}

// STL to DataFrame conversion functions
DataFrame point_rows_to_dataframe(const std::vector<PointRow>& rows) {
  int n = rows.size();
  
  if (n == 0) {
    return DataFrame::create(
      Named("contig") = CharacterVector(),
      Named("coord") = IntegerVector(),
      Named("vcoord") = NumericVector(),
      Named("input_index") = IntegerVector(),
      Named("stringsAsFactors") = false
    );
  }
  
  CharacterVector contigs(n);
  IntegerVector coords(n);
  NumericVector vcoords(n);
  IntegerVector input_indices(n);
  
  for (int i = 0; i < n; i++) {
    contigs[i] = rows[i].contig;
    coords[i] = static_cast<int>(rows[i].coord);
    vcoords[i] = rows[i].vcoord;
    input_indices[i] = rows[i].input_index;
  }
  
  return DataFrame::create(
    Named("contig") = contigs,
    Named("coord") = coords,
    Named("vcoord") = vcoords,
    Named("input_index") = input_indices,
    Named("stringsAsFactors") = false
  );
}

DataFrame interval_rows_to_dataframe(const std::vector<IntervalRow>& rows) {
  int n = rows.size();
  
  if (n == 0) {
    return DataFrame::create(
      Named("contig") = CharacterVector(),
      Named("start") = IntegerVector(),
      Named("end") = IntegerVector(),
      Named("vstart") = NumericVector(),
      Named("vend") = NumericVector(),
      Named("trim_left") = LogicalVector(),
      Named("trim_right") = LogicalVector(),
      Named("input_index") = IntegerVector(),
      Named("stringsAsFactors") = false
    );
  }
  
  CharacterVector contigs(n);
  IntegerVector starts(n);
  IntegerVector ends(n);
  NumericVector vstarts(n);
  NumericVector vends(n);
  LogicalVector trim_left(n);
  LogicalVector trim_right(n);
  IntegerVector input_indices(n);
  
  for (int i = 0; i < n; i++) {
    contigs[i] = rows[i].contig;
    starts[i] = static_cast<int>(rows[i].start);
    ends[i] = static_cast<int>(rows[i].end);
    vstarts[i] = rows[i].vstart;
    vends[i] = rows[i].vend;
    trim_left[i] = rows[i].trim_left;
    trim_right[i] = rows[i].trim_right;
    input_indices[i] = rows[i].input_index;
  }
  
  return DataFrame::create(
    Named("contig") = contigs,
    Named("start") = starts,
    Named("end") = ends,
    Named("vstart") = vstarts,
    Named("vend") = vends,
    Named("trim_left") = trim_left,
    Named("trim_right") = trim_right,
    Named("input_index") = input_indices,
    Named("stringsAsFactors") = false
  );
}

DataFrame plotted_segments_to_dataframe(const std::vector<PlottedSegment>& segments) {
  int n = segments.size();
  
  if (n == 0) {
    // Return a proper empty DataFrame with the expected columns
    return DataFrame::create(
      Named("segment") = CharacterVector(),
      Named("contig") = CharacterVector(),
      Named("start") = NumericVector(),
      Named("end") = NumericVector(),
      Named("vstart") = NumericVector(),
      Named("vend") = NumericVector(),
      Named("contig_start") = IntegerVector(),
      Named("contig_end") = IntegerVector(),
      Named("stringsAsFactors") = false
    );
  }
  
  CharacterVector segment_ids(n);
  CharacterVector contigs(n);
  NumericVector vstarts(n);
  NumericVector vends(n);
  IntegerVector contig_starts(n);
  IntegerVector contig_ends(n);
  
  for (int i = 0; i < n; i++) {
    segment_ids[i] = segments[i].segment_id;
    contigs[i] = segments[i].contig;
    vstarts[i] = static_cast<double>(segments[i].vcoord_start);
    vends[i] = static_cast<double>(segments[i].vcoord_end);
    contig_starts[i] = static_cast<int>(segments[i].contig_start);
    contig_ends[i] = static_cast<int>(segments[i].contig_end);
  }
  
  return DataFrame::create(
    Named("segment") = segment_ids,
    Named("contig") = contigs,
    Named("start") = vstarts,  // vcoord start for compatibility
    Named("end") = vends,      // vcoord end for compatibility
    Named("vstart") = vstarts,
    Named("vend") = vends,
    Named("contig_start") = contig_starts,
    Named("contig_end") = contig_ends,
    Named("stringsAsFactors") = false
  );
}

// For get_view_intervals: vstart, vend, contig, start, end, segment_ids, n_segments
DataFrame view_intervals_to_dataframe(const std::vector<IntervalRow>& rows) {
  int n = rows.size();
  
  if (n == 0) {
    return DataFrame::create(
      Named("vstart") = NumericVector(),
      Named("vend") = NumericVector(),
      Named("contig") = CharacterVector(),
      Named("start") = IntegerVector(),
      Named("end") = IntegerVector(),
      Named("segment_ids") = CharacterVector(),
      Named("n_segments") = IntegerVector(),
      Named("stringsAsFactors") = false
    );
  }
  
  NumericVector vstarts(n);
  NumericVector vends(n);
  CharacterVector contigs(n);
  IntegerVector starts(n);
  IntegerVector ends(n);
  CharacterVector segment_ids(n);
  IntegerVector n_segments(n);
  
  for (int i = 0; i < n; i++) {
    vstarts[i] = rows[i].vstart;
    vends[i] = rows[i].vend;
    contigs[i] = rows[i].contig;
    starts[i] = static_cast<int>(rows[i].start);
    ends[i] = static_cast<int>(rows[i].end);
    segment_ids[i] = rows[i].segment_ids;
    n_segments[i] = rows[i].n_segments;
  }
  
  return DataFrame::create(
    Named("vstart") = vstarts,
    Named("vend") = vends,
    Named("contig") = contigs,
    Named("start") = starts,
    Named("end") = ends,
    Named("segment_ids") = segment_ids,
    Named("n_segments") = n_segments,
    Named("stringsAsFactors") = false
  );
}
}

// [[Rcpp::export]]
ContextPtr context_create() {
  Context* ctx = new Context();  // Empty constructor
  return ContextPtr(ctx, true);
}

// [[Rcpp::export]]
void context_set_tables(ContextPtr ctx, DataFrame segment_table, DataFrame contig_table) {
  try {
    Context* ptr = require_context(ctx, "context_set_tables");
    std::vector<SegmentRow> segment_rows = dataframe_to_segment_rows(segment_table);
    std::vector<ContigRow> contig_rows = dataframe_to_contig_rows(contig_table);
    ptr->set_tables(segment_rows, contig_rows);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_set_tables: %s", e.what());
  }
}

// [[Rcpp::export]]
void context_update_selected_segments(ContextPtr ctx, DataFrame selected_segments, 
                                      NumericVector zoom = NumericVector()) {
  try {
    Context* ptr = require_context(ctx, "context_update_selected_segments");
    std::vector<SegmentRow> segment_rows = dataframe_to_segment_rows(selected_segments);
    std::vector<double> zoom_vec = to_std_vector(zoom);
    ptr->update_selected_segments(segment_rows, zoom_vec);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_update_selected_segments: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_contig2view_point(ContextPtr ctx, DataFrame df) {
  try {
    Context* ptr = require_context(ctx, "context_contig2view_point");
    std::vector<PointRow> input = dataframe_to_point_rows(df, false);
    std::vector<PointRow> result = ptr->contig2view_point(input);
    return point_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_contig2view_point: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_view2contig_point(ContextPtr ctx, DataFrame df) {
  try {
    Context* ptr = require_context(ctx, "context_view2contig_point");
    std::vector<PointRow> input = dataframe_to_point_rows(df, true);
    std::vector<PointRow> result = ptr->view2contig_point(input);
    return point_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_view2contig_point: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_contig2view_interval(ContextPtr ctx, DataFrame df, 
                                        String crossing_intervals = "drop") {
  try {
    Context* ptr = require_context(ctx, "context_contig2view_interval");
    std::vector<IntervalRow> input = dataframe_to_interval_rows(df, false);
    std::vector<IntervalRow> result = ptr->contig2view_interval(input, std::string(crossing_intervals));
    return interval_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_contig2view_interval: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_view2contig_interval(ContextPtr ctx, DataFrame df) {
  try {
    Context* ptr = require_context(ctx, "context_view2contig_interval");
    std::vector<IntervalRow> input = dataframe_to_interval_rows(df, true);
    std::vector<IntervalRow> result = ptr->view2contig_interval(input);
    return interval_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_view2contig_interval: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_get_view_intervals(ContextPtr ctx, bool limit_to_zoom = true, bool merge_adjacent = true) {
  try {
    Context* ptr = require_context(ctx, "context_get_view_intervals");
    std::vector<IntervalRow> result = ptr->get_view_intervals(limit_to_zoom, merge_adjacent);
    return view_intervals_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_get_view_intervals: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_filter_coords(ContextPtr ctx, DataFrame df, 
                                NumericVector xlim = NumericVector()) {
  try {
    Context* ptr = require_context(ctx, "context_filter_coords");
    std::vector<PointRow> input = dataframe_to_point_rows(df, false);
    std::vector<double> xlim_vec = to_std_vector(xlim);
    std::vector<PointRow> result = ptr->filter_coords(input, xlim_vec);
    return point_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_filter_coords: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_filter_intervals(ContextPtr ctx, DataFrame df, 
                                    NumericVector xlim = NumericVector(),
                                    bool merge_adjacent = false) {
  try {
    Context* ptr = require_context(ctx, "context_filter_intervals");
    std::vector<IntervalRow> input = dataframe_to_interval_rows(df, false);
    std::vector<double> xlim_vec = to_std_vector(xlim);
    std::vector<IntervalRow> result = ptr->filter_intervals(input, xlim_vec, merge_adjacent);
    return interval_rows_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_filter_intervals: %s", e.what());
  }
}

// [[Rcpp::export]]
LogicalVector context_coords_in_view(ContextPtr ctx, DataFrame df, bool limit_to_zoom = false) {
  try {
    Context* ptr = require_context(ctx, "context_coords_in_view");
    std::vector<PointRow> input = dataframe_to_point_rows(df, false);
    std::vector<bool> result = ptr->coords_in_view(input, limit_to_zoom);
    return Rcpp::wrap(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_coords_in_view: %s", e.what());
  }
}

// [[Rcpp::export]]
DataFrame context_get_plotted_segments(ContextPtr ctx) {
  try {
    Context* ptr = require_context(ctx, "context_get_plotted_segments");
    std::vector<PlottedSegment> result = ptr->get_plotted_segments();
    return plotted_segments_to_dataframe(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_get_plotted_segments: %s", e.what());
  }
}

// [[Rcpp::export]]
NumericVector context_get_xlim(ContextPtr ctx) {
  try {
    Context* ptr = require_context(ctx, "context_get_xlim");
    std::vector<double> xlim = ptr->get_xlim();
    NumericVector result(2);
    if (xlim.size() == 2) {
      result[0] = xlim[0];
      result[1] = xlim[1];
    } else {
      result[0] = NA_REAL;
      result[1] = NA_REAL;
    }
    return result;
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_get_xlim: %s", e.what());
  }
}

// [[Rcpp::export]]
CharacterVector context_get_segment_ids(ContextPtr ctx, bool limit_to_zoom = false) {
  try {
    Context* ptr = require_context(ctx, "context_get_segment_ids");
    std::vector<std::string> result = ptr->get_segment_ids(limit_to_zoom);
    return wrap(result);
  } catch (const std::exception& e) {
    Rcpp::stop("Error in context_get_segment_ids: %s", e.what());
  }
}

