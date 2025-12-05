#include "Context.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <stdexcept>
#include <sstream>

Context::Context() {
  // Do nothing - tables will be set via set_tables()
}

void Context::set_tables(const std::vector<SegmentRow>& segment_table, 
                         const std::vector<ContigRow>& contig_table) {
  build_full_segments(segment_table);
  build_contigs(contig_table);
}

void Context::build_full_segments(const std::vector<SegmentRow>& segment_table) {
  size_t n = segment_table.size();
  full_segments_.clear();
  full_segments_.reserve(n);
  segment_id_to_index_.clear();
  
  int64_t acoord_cumsum = 0;
  
  for (size_t i = 0; i < n; i++) {
    const SegmentRow& row = segment_table[i];
    Segment seg;
    seg.segment_id = row.segment_id;
    seg.contig = row.contig;
    seg.start = row.start;
    seg.end = row.end;
    seg.length = seg.end - seg.start + 1;
    seg.acoord_start = acoord_cumsum;
    
    full_segments_.push_back(seg);
    segment_id_to_index_[seg.segment_id] = i;
    
    acoord_cumsum += seg.length;
  }
  
  // Build prefix sum for binary search (exclusive ends)
  all_segment_cs.resize(n);
  for (size_t i = 0; i < n; i++) {
    all_segment_cs[i] = full_segments_[i].acoord_start + full_segments_[i].length;
  }
}

void Context::build_contigs(const std::vector<ContigRow>& contig_table) {
  size_t n = contig_table.size();
  contigs_.clear();
  contigs_.reserve(n);
  contig_id_to_index_.clear();
  
  int64_t acoord_cumsum = 0;
  
  // Build contig table sequentially in acoord space (independent of segments)
  for (size_t i = 0; i < n; i++) {
    const ContigRow& row = contig_table[i];
    Contig ctg;
    ctg.contig_id = row.contig_id;
    ctg.length = row.length;
    ctg.acoord_start = acoord_cumsum;
    
    contigs_.push_back(ctg);
    contig_id_to_index_[ctg.contig_id] = i;
    
    acoord_cumsum += ctg.length;
  }
  
  // Build prefix sum for binary search (exclusive ends)
  all_contig_cs.resize(n);
  for (size_t i = 0; i < n; i++) {
    all_contig_cs[i] = contigs_[i].acoord_start + contigs_[i].length;
  }
}

void Context::update_selected_segments(const std::vector<SegmentRow>& selected_segments, 
                                       const std::vector<double>& zoom) {
  update_plotted_segments(selected_segments);
  update_xlim(zoom);
}

void Context::update_plotted_segments(const std::vector<SegmentRow>& selected_segments) {
  size_t n = selected_segments.size();
  plotted_segments_.clear();
  plotted_segments_.reserve(n);
  plotted_segment_id_to_index_.clear();
  
  int64_t vcoord_cumsum = 0;
  
  for (size_t i = 0; i < n; i++) {
    const SegmentRow& row = selected_segments[i];
    std::string seg_id = row.segment_id;
    auto it = segment_id_to_index_.find(seg_id);
    if (it == segment_id_to_index_.end()) {
      continue; // skip invalid segments
    }
    
    PlottedSegment pseg;
    pseg.segment_id = seg_id;
    pseg.contig = row.contig;
    pseg.contig_start = row.start;
    pseg.contig_end = row.end;
    pseg.strand = row.strand;
    
    int64_t length = pseg.contig_end - pseg.contig_start + 1;
    pseg.vcoord_start = vcoord_cumsum;
    pseg.vcoord_end = vcoord_cumsum + length - 1;
    
    plotted_segments_.push_back(pseg);
    plotted_segment_id_to_index_[seg_id] = plotted_segments_.size() - 1;
    vcoord_cumsum += length;
  }
  
  // Build prefix sum for binary search (exclusive ends)
  plotted_segment_cs.resize(plotted_segments_.size());
  for (size_t i = 0; i < plotted_segments_.size(); i++) {
    plotted_segment_cs[i] = plotted_segments_[i].vcoord_end + 1;
  }
}

void Context::update_xlim(const std::vector<double>& zoom) {
  if (plotted_segments_.empty()) {
    xlim_ = std::vector<double>(2, 0.0);
    return;
  }
  
  if (zoom.size() == 2) {
    xlim_ = zoom;
  } else {
    xlim_.resize(2);
    xlim_[0] = static_cast<double>(plotted_segments_[0].vcoord_start);
    xlim_[1] = static_cast<double>(plotted_segments_.back().vcoord_end);
  }
}

int64_t Context::contig2acoord(const std::string& contig, int64_t coord) {
  auto it = contig_id_to_index_.find(contig);
  if (it == contig_id_to_index_.end()) {
    std::ostringstream oss;
    oss << "error: contig not found in context: '" << contig << "'";
    throw std::runtime_error(oss.str());
  }
  
  const Contig& ctg = contigs_[it->second];
  // coord is 1-based on contig, convert to 0-based acoord
  return ctg.acoord_start + (coord - 1);
}

std::pair<std::string, int64_t> Context::acoord2contig(int64_t acoord) {
  if (all_contig_cs.empty() || acoord < 0) {
    std::ostringstream oss;
    oss << "error: acoord out of range: " << acoord;
    throw std::runtime_error(oss.str());
  }
  
  // Binary search on contigs: find first contig prefix > acoord
  auto it = std::upper_bound(all_contig_cs.begin(), all_contig_cs.end(), acoord);
  size_t contig_idx = it - all_contig_cs.begin();
  
  if (contig_idx >= contigs_.size()) {
    std::ostringstream oss;
    oss << "error: acoord " << acoord << " not found in any contig";
    if (!all_contig_cs.empty()) {
      oss << " (valid range: [0, " << all_contig_cs.back() << "))";
    }
    throw std::runtime_error(oss.str());
  }
  
  const Contig& ctg = contigs_[contig_idx];
  
  // Direct conversion: acoord is 0-based, convert to 1-based contig coordinate
  int64_t local_acoord = acoord - ctg.acoord_start; // 0-based within contig
  int64_t contig_coord = local_acoord + 1; // convert to 1-based
  
  return std::make_pair(ctg.contig_id, contig_coord);
}

std::pair<std::string, int64_t> Context::acoord2segment(int64_t acoord) {
  if (all_segment_cs.empty() || acoord < 0 || acoord >= all_segment_cs.back()) {
    std::ostringstream oss;
    oss << "error: acoord out of range: " << acoord;
    throw std::runtime_error(oss.str());
  }
  
  // Binary search: find first prefix > acoord
  auto it = std::upper_bound(all_segment_cs.begin(), all_segment_cs.end(), acoord);
  size_t idx = it - all_segment_cs.begin();
  
  const Segment& seg = full_segments_[idx];
  int64_t local_acoord = acoord - seg.acoord_start; // 0-based within segment
  int64_t coord = local_acoord + 1; // convert to 1-based
  
  return std::make_pair(seg.segment_id, coord);
}

int64_t Context::segment2acoord(const std::string& segment, int64_t coord) {
  auto it = segment_id_to_index_.find(segment);
  if (it == segment_id_to_index_.end()) {
    std::ostringstream oss;
    oss << "error: segment not found in context: '" << segment << "'";
    throw std::runtime_error(oss.str());
  }
  return full_segments_[it->second].acoord_start + coord - 1; // coord is 1-based
}

int64_t Context::segment2vcoord(const std::string& segment, int64_t coord) {
  // coord is 1-based within segment space
  auto it = plotted_segment_id_to_index_.find(segment);
  if (it == plotted_segment_id_to_index_.end()) {
    std::ostringstream oss;
    oss << "error: segment not in plotted segments: '" << segment << "'";
    throw std::runtime_error(oss.str());
  }
  
  const PlottedSegment& pseg = plotted_segments_[it->second];
  int64_t local_coord = coord - 1; // convert to 0-based
  return segment_local_to_vcoord(pseg, local_coord);
}

std::pair<std::string, int64_t> Context::vcoord2segment(int64_t vcoord) {
  if (plotted_segment_cs.empty() || vcoord < 0 || vcoord >= plotted_segment_cs.back()) {
    std::ostringstream oss;
    oss << "error: vcoord out of range: " << vcoord << " (valid range depends on plotted segments)";
    throw std::runtime_error(oss.str());
  }
  
  // Binary search: find first prefix > vcoord
  auto it = std::upper_bound(plotted_segment_cs.begin(), plotted_segment_cs.end(), vcoord);
  size_t idx = it - plotted_segment_cs.begin();
  
  const PlottedSegment& pseg = plotted_segments_[idx];
  int64_t local_coord = vcoord_to_segment_local(pseg, vcoord);
  int64_t coord = local_coord + 1; // convert to 1-based
  
  return std::make_pair(pseg.segment_id, coord);
}

int64_t Context::segment_local_to_vcoord(const PlottedSegment& pseg, int64_t local_coord) {
  // local_coord is 0-based offset within segment (0 to length-1)
  // for '+' strand: vcoord = vcoord_start + local_coord
  // for '-' strand: vcoord = vcoord_end - local_coord (reversed)
  if (pseg.strand == '-') {
    return pseg.vcoord_end - local_coord;
  }
  return pseg.vcoord_start + local_coord;
}

int64_t Context::vcoord_to_segment_local(const PlottedSegment& pseg, int64_t vcoord) {
  // reverse of segment_local_to_vcoord
  // returns 0-based local coordinate within segment
  if (pseg.strand == '-') {
    return pseg.vcoord_end - vcoord;
  }
  return vcoord - pseg.vcoord_start;
}

std::vector<PointRow> Context::contig2view_point(const std::vector<PointRow>& input) {
  std::vector<PointRow> result;
  result.reserve(input.size());
  
  for (size_t i = 0; i < input.size(); i++) {
    const PointRow& row = input[i];
    PointRow out_row = row;
    out_row.input_index = static_cast<int>(i + 1);  // 1-based
    
    // contig -> acoord -> segment -> vcoord
    int64_t acoord = contig2acoord(row.contig, row.coord);
    auto seg_info = acoord2segment(acoord);
    int64_t vcoord = segment2vcoord(seg_info.first, seg_info.second);
    
    out_row.vcoord = static_cast<double>(vcoord);
    result.push_back(out_row);
  }
  
  return result;
}

std::vector<PointRow> Context::view2contig_point(const std::vector<PointRow>& input) {
  std::vector<PointRow> result;
  result.reserve(input.size());
  
  for (size_t i = 0; i < input.size(); i++) {
    const PointRow& row = input[i];
    PointRow out_row = row;
    out_row.input_index = static_cast<int>(i + 1);  // 1-based
    int64_t vcoord = static_cast<int64_t>(row.vcoord);
    
    // vcoord -> segment -> acoord -> contig
    auto seg_info = vcoord2segment(vcoord);
    int64_t acoord = segment2acoord(seg_info.first, seg_info.second);
    auto contig_info = acoord2contig(acoord);
    
    out_row.contig = contig_info.first;
    out_row.coord = contig_info.second;
    result.push_back(out_row);
  }
  
  return result;
}

std::vector<IntervalRow> Context::contig2view_interval(const std::vector<IntervalRow>& input, 
                                                        const std::string& crossing_intervals) {
  std::vector<IntervalRow> result;
  
  for (size_t i = 0; i < input.size(); i++) {
    const IntervalRow& row = input[i];
    int input_idx = static_cast<int>(i + 1);  // 1-based
    std::string contig = row.contig;
    int64_t start = row.start;
    int64_t end = row.end;
    
    // Find all plotted segments on this contig that overlap the interval
    std::vector<const PlottedSegment*> overlapping_segments;
    for (const auto& pseg : plotted_segments_) {
      if (pseg.contig == contig) {
        // Check if segment overlaps with interval [start, end]
        if (pseg.contig_end >= start && pseg.contig_start <= end) {
          overlapping_segments.push_back(&pseg);
        }
      }
    }
    
    if (overlapping_segments.empty()) {
      continue; // No overlapping segments, skip
    }
    
    if (crossing_intervals == "split") {
      // Return one interval per overlapping segment
      for (const auto* pseg : overlapping_segments) {
        // Compute intersection on contig
        int64_t intersect_start = std::max(start, pseg->contig_start);
        int64_t intersect_end = std::min(end, pseg->contig_end);
        
        // Transform intersection to vcoord using strand-aware helper
        int64_t local_start = intersect_start - pseg->contig_start;
        int64_t local_end = intersect_end - pseg->contig_start;
        int64_t v1 = segment_local_to_vcoord(*pseg, local_start);
        int64_t v2 = segment_local_to_vcoord(*pseg, local_end);
        int64_t vstart = std::min(v1, v2);
        int64_t vend = std::max(v1, v2);
        
        // Check for trimming (swap for reversed strand)
        bool trim_contig_left = (intersect_start > start);
        bool trim_contig_right = (intersect_end < end);
        bool trim_l = (pseg->strand == '-') ? trim_contig_right : trim_contig_left;
        bool trim_r = (pseg->strand == '-') ? trim_contig_left : trim_contig_right;
        
        IntervalRow out_row;
        out_row.contig = contig;
        out_row.start = row.start;
        out_row.end = row.end;
        out_row.vstart = static_cast<double>(vstart);
        out_row.vend = static_cast<double>(vend);
        out_row.trim_left = trim_l;
        out_row.trim_right = trim_r;
        out_row.n_segments = 0;
        out_row.input_index = input_idx;
        out_row.segment_strand = pseg->strand;
        result.push_back(out_row);
      }
    } else {
      // "drop" mode: only process if interval is within a single segment
      if (overlapping_segments.size() == 1) {
        const PlottedSegment* pseg = overlapping_segments[0];
        
        // Check if entire interval is within this segment
        if (start >= pseg->contig_start && end <= pseg->contig_end) {
          // Transform to vcoord using strand-aware helper
          int64_t local_start = start - pseg->contig_start;
          int64_t local_end = end - pseg->contig_start;
          int64_t v1 = segment_local_to_vcoord(*pseg, local_start);
          int64_t v2 = segment_local_to_vcoord(*pseg, local_end);
          int64_t vstart = std::min(v1, v2);
          int64_t vend = std::max(v1, v2);
          
          // No trimming since interval is fully within segment
          IntervalRow out_row;
          out_row.contig = contig;
          out_row.start = row.start;
          out_row.end = row.end;
          out_row.vstart = static_cast<double>(vstart);
          out_row.vend = static_cast<double>(vend);
          out_row.trim_left = false;
          out_row.trim_right = false;
          out_row.n_segments = 0;
          out_row.input_index = input_idx;
          out_row.segment_strand = pseg->strand;
          result.push_back(out_row);
        } else {
          // Interval extends beyond segment boundaries, check trimming
          int64_t intersect_start = std::max(start, pseg->contig_start);
          int64_t intersect_end = std::min(end, pseg->contig_end);
          
          int64_t local_start = intersect_start - pseg->contig_start;
          int64_t local_end = intersect_end - pseg->contig_start;
          int64_t v1 = segment_local_to_vcoord(*pseg, local_start);
          int64_t v2 = segment_local_to_vcoord(*pseg, local_end);
          int64_t vstart = std::min(v1, v2);
          int64_t vend = std::max(v1, v2);
          
          bool trim_contig_left = (intersect_start > start);
          bool trim_contig_right = (intersect_end < end);
          bool trim_l = (pseg->strand == '-') ? trim_contig_right : trim_contig_left;
          bool trim_r = (pseg->strand == '-') ? trim_contig_left : trim_contig_right;
          
          IntervalRow out_row;
          out_row.contig = contig;
          out_row.start = row.start;
          out_row.end = row.end;
          out_row.vstart = static_cast<double>(vstart);
          out_row.vend = static_cast<double>(vend);
          out_row.trim_left = trim_l;
          out_row.trim_right = trim_r;
          out_row.n_segments = 0;
          out_row.input_index = input_idx;
          out_row.segment_strand = pseg->strand;
          result.push_back(out_row);
        }
      }
      // else: interval crosses multiple segments, drop it
    }
  }
  
  return result;
}

std::vector<IntervalRow> Context::view2contig_interval(const std::vector<IntervalRow>& input) {
  std::vector<IntervalRow> result;
  result.reserve(input.size());
  
  for (size_t i = 0; i < input.size(); i++) {
    const IntervalRow& row = input[i];
    
    IntervalRow out_row = row;
    out_row.input_index = static_cast<int>(i + 1);  // 1-based
    int64_t vstart = static_cast<int64_t>(row.vstart);
    int64_t vend = static_cast<int64_t>(row.vend);
    
    // vcoord -> segment -> acoord -> contig
    auto seg_start_info = vcoord2segment(vstart);
    
    auto seg_end_info = vcoord2segment(vend);
    
    int64_t acoord_start = segment2acoord(seg_start_info.first, seg_start_info.second);
    int64_t acoord_end = segment2acoord(seg_end_info.first, seg_end_info.second);
    
    auto contig_start_info = acoord2contig(acoord_start);
    auto contig_end_info = acoord2contig(acoord_end);
    
    out_row.contig = contig_start_info.first;
    out_row.start = contig_start_info.second;
    out_row.end = contig_end_info.second;
    
    result.push_back(out_row);
  }
  
  return result;
}

std::vector<IntervalRow> Context::get_view_intervals(bool limit_to_zoom, bool merge_adjacent) {
  std::vector<IntervalRow> intervals;
  
  if (plotted_segments_.empty()) {
    return intervals;
  }
  
  // Step 1: Build intervals from plotted segments, optionally clipped to xlim
  for (const auto& pseg : plotted_segments_) {
    int64_t vstart = pseg.vcoord_start;
    int64_t vend = pseg.vcoord_end;
    int64_t cstart = pseg.contig_start;
    int64_t cend = pseg.contig_end;
    
    if (limit_to_zoom && xlim_.size() == 2) {
      int64_t xlim_start = static_cast<int64_t>(xlim_[0]);
      int64_t xlim_end = static_cast<int64_t>(xlim_[1]);
      
      // Check overlap
      if (vend < xlim_start || vstart > xlim_end) {
        continue;  // No overlap, skip
      }
      
      // Clip vcoords to xlim
      int64_t new_vstart = std::max(vstart, xlim_start);
      int64_t new_vend = std::min(vend, xlim_end);
      
      // Adjust contig coords based on strand
      int64_t offset_start = new_vstart - vstart;
      int64_t offset_end = vend - new_vend;
      int64_t new_cstart, new_cend;
      if (pseg.strand == '-') {
        // reversed: vcoord increases as contig coord decreases
        new_cstart = cstart + offset_end;
        new_cend = cend - offset_start;
      } else {
        new_cstart = cstart + offset_start;
        new_cend = cend - offset_end;
      }
      
      vstart = new_vstart;
      vend = new_vend;
      cstart = new_cstart;
      cend = new_cend;
    }
    
    IntervalRow row;
    row.vstart = static_cast<double>(vstart);
    row.vend = static_cast<double>(vend);
    row.contig = pseg.contig;
    row.start = cstart;
    row.end = cend;
    row.segment_ids = pseg.segment_id;
    row.n_segments = 1;
    row.trim_left = false;
    row.trim_right = false;
    row.segment_strand = pseg.strand;
    intervals.push_back(row);
  }
  
  if (!merge_adjacent || intervals.empty()) {
    return intervals;
  }
  
  // Step 2: Merge adjacent intervals on same contig with same strand
  
  std::vector<IntervalRow> merged;
  IntervalRow current = intervals[0];
  
  for (size_t i = 1; i < intervals.size(); i++) {
    const IntervalRow& next = intervals[i];
    
    // Check if same contig and same strand
    bool same_contig = (next.contig == current.contig);
    bool same_strand = (next.segment_strand == current.segment_strand);
    
    // Check adjacency based on strand
    bool adjacent = false;
    if (same_contig && same_strand) {
      if (current.segment_strand == '+') {
        // plus strand: segments in contig order
        adjacent = (next.start == current.end + 1);
      } else {
        // minus strand: segments in reverse contig order (current has higher coords)
        adjacent = (current.start == next.end + 1);
      }
    }
    
    if (adjacent) {
      // Merge: extend current interval
      if (current.segment_strand == '+') {
        current.end = next.end;
      } else {
        current.start = next.start;
      }
      current.vend = next.vend;
      current.segment_ids += "," + next.segment_ids;
      current.n_segments += next.n_segments;
    } else {
      // Not adjacent, push current and start new
      merged.push_back(current);
      current = next;
    }
  }
  // Push last interval
  merged.push_back(current);
  
  return merged;
}

std::vector<PointRow> Context::filter_coords(const std::vector<PointRow>& input, 
                                             const std::vector<double>& xlim) {
  if (input.empty()) {
    return std::vector<PointRow>();
  }
  
  std::vector<PointRow> result;
  result.reserve(input.size());
  
  for (size_t i = 0; i < input.size(); i++) {
    const PointRow& row = input[i];
    
    // Convert contig coord to acoord
    auto contig_it = contig_id_to_index_.find(row.contig);
    if (contig_it == contig_id_to_index_.end()) {
      continue; // skip unknown contigs
    }
    
    int64_t acoord = contig2acoord(row.contig, row.coord);
    
    // Find which segment contains this acoord
    if (all_segment_cs.empty() || acoord < 0 || acoord >= all_segment_cs.back()) {
      continue; // skip out of range
    }
    auto seg_it = std::upper_bound(all_segment_cs.begin(), all_segment_cs.end(), acoord);
    size_t seg_idx = seg_it - all_segment_cs.begin();
    const Segment& seg = full_segments_[seg_idx];
    
    // Check if this segment is plotted
    auto plotted_it = plotted_segment_id_to_index_.find(seg.segment_id);
    if (plotted_it == plotted_segment_id_to_index_.end()) {
      continue; // skip if segment not plotted
    }
    
    // Convert to vcoord using strand-aware helper
    const PlottedSegment& pseg = plotted_segments_[plotted_it->second];
    int64_t segment_local = acoord - seg.acoord_start; // 0-based within segment
    int64_t vcoord = segment_local_to_vcoord(pseg, segment_local);
    
    // Check xlim
    if (xlim.size() == 2) {
      if (vcoord < static_cast<int64_t>(xlim[0]) || vcoord > static_cast<int64_t>(xlim[1])) {
        continue;
      }
    }
    
    PointRow out_row = row;
    out_row.input_index = static_cast<int>(i + 1);  // 1-based
    out_row.vcoord = static_cast<double>(vcoord);
    result.push_back(out_row);
  }
  
  return result;
}

std::vector<bool> Context::coords_in_view(const std::vector<PointRow>& input, bool limit_to_zoom) {
  std::vector<bool> result(input.size(), false);
  
  for (size_t i = 0; i < input.size(); i++) {
    const PointRow& row = input[i];
    
    auto contig_it = contig_id_to_index_.find(row.contig);
    if (contig_it == contig_id_to_index_.end()) {
      continue;
    }
    
    int64_t acoord = contig2acoord(row.contig, row.coord);
    
    if (all_segment_cs.empty() || acoord < 0 || acoord >= all_segment_cs.back()) {
      continue;
    }
    
    auto seg_it = std::upper_bound(all_segment_cs.begin(), all_segment_cs.end(), acoord);
    size_t seg_idx = seg_it - all_segment_cs.begin();
    const Segment& seg = full_segments_[seg_idx];
    
    auto plotted_it = plotted_segment_id_to_index_.find(seg.segment_id);
    if (plotted_it == plotted_segment_id_to_index_.end()) {
      continue;
    }
    
    if (limit_to_zoom && xlim_.size() == 2) {
      const PlottedSegment& pseg = plotted_segments_[plotted_it->second];
      int64_t segment_local = acoord - seg.acoord_start;
      int64_t vcoord = segment_local_to_vcoord(pseg, segment_local);
      if (vcoord < static_cast<int64_t>(xlim_[0]) || vcoord > static_cast<int64_t>(xlim_[1])) {
        continue;
      }
    }
    
    result[i] = true;
  }
  
  return result;
}

std::vector<IntervalRow> Context::filter_intervals(const std::vector<IntervalRow>& input, 
                                                    const std::vector<double>& xlim,
                                                    bool merge_adjacent) {
  if (input.empty()) {
    return std::vector<IntervalRow>();
  }
  
  // Transform intervals to vcoord (split across segments)
  std::vector<IntervalRow> with_vcoord = contig2view_interval(input, "split");
  
  if (with_vcoord.empty()) {
    return std::vector<IntervalRow>();
  }
  
  // Filter by xlim if provided
  std::vector<IntervalRow> filtered;
  filtered.reserve(with_vcoord.size());
  
  // Also filter to only contigs in plotted segments
  std::unordered_set<std::string> valid_contigs;
  for (const auto& pseg : plotted_segments_) {
    valid_contigs.insert(pseg.contig);
  }
  
  for (const auto& row : with_vcoord) {
    bool keep = true;
    
    // Check xlim overlap
    if (xlim.size() == 2) {
      if (row.vend < xlim[0] || row.vstart > xlim[1]) {
        keep = false;
      }
    }
    
    // Check valid contigs
    if (keep && valid_contigs.find(row.contig) == valid_contigs.end()) {
      keep = false;
    }
    
    if (keep) {
      filtered.push_back(row);
    }
  }
  
  if (!merge_adjacent || filtered.empty()) {
    return filtered;
  }
  
  // Sort by input_index first, then by vstart, so intervals from same input row are grouped
  std::sort(filtered.begin(), filtered.end(), [](const IntervalRow& a, const IntervalRow& b) {
    if (a.input_index != b.input_index) {
      return a.input_index < b.input_index;
    }
    return a.vstart < b.vstart;
  });
  
  // Merge adjacent intervals (adjacent in vcoord space, same strand, same input_index)
  std::vector<IntervalRow> merged;
  IntervalRow current = filtered[0];
  
  for (size_t i = 1; i < filtered.size(); i++) {
    const IntervalRow& next = filtered[i];
    
    // Check if adjacent in vcoord space, same strand, and same input_index
    int64_t current_vend = static_cast<int64_t>(current.vend);
    int64_t next_vstart = static_cast<int64_t>(next.vstart);
    bool same_strand = (next.segment_strand == current.segment_strand);
    bool same_input = (next.input_index == current.input_index);
    
    if (same_input && same_strand && next_vstart == current_vend + 1) {
      current.end = next.end;
      current.vend = next.vend;
      current.trim_right = next.trim_right;
      current.n_segments += next.n_segments;
    } else {
      merged.push_back(current);
      current = next;
    }
  }
  merged.push_back(current);
  
  return merged;
}

std::vector<PlottedSegment> Context::get_plotted_segments() {
  return plotted_segments_;
}

std::vector<double> Context::get_xlim() {
  return xlim_;
}

std::vector<std::string> Context::get_segment_ids(bool limit_to_zoom) {
  std::vector<std::string> result;
  
  if (!limit_to_zoom) {
    // return all plotted segment IDs
    for (const auto& seg : plotted_segments_) {
      result.push_back(seg.segment_id);
    }
    return result;
  }
  
  // limit to segments visible in current zoom
  if (xlim_.size() != 2) {
    // no zoom set, return all
    for (const auto& seg : plotted_segments_) {
      result.push_back(seg.segment_id);
    }
    return result;
  }
  
  double zoom_start = xlim_[0];
  double zoom_end = xlim_[1];
  
  for (const auto& seg : plotted_segments_) {
    // check if segment overlaps with zoom range
    if (seg.vcoord_end >= zoom_start && seg.vcoord_start <= zoom_end) {
      result.push_back(seg.segment_id);
    }
  }
  
  return result;
}


