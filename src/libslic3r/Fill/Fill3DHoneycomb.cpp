///|/ Copyright (c) Prusa Research 2016 - 2021 Vojtěch Bubník @bubnikv
///|/ Copyright (c) SuperSlicer 2019 Remi Durand @supermerill
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include <algorithm>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstddef>

#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "Fill3DHoneycomb.hpp"
#include "libslic3r/BoundingBox.hpp"
#include "libslic3r/Fill/FillBase.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/libslic3r.h"

namespace Slic3r {

// sign function
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
  
/*
Creates a contiguous sequence of points at a specified height that make
up a horizontal slice of the edges of a space filling truncated
octahedron tesselation. The octahedrons are oriented so that the
square faces are in the horizontal plane with edges parallel to the X
and Y axes.

Credits: David Eccles (gringer).
*/

// triangular wave function
// this has period (gridSize * 2), and amplitude (gridSize / 2),
// with triWave(pos = 0) = 0
static coordf_t triWave(coordf_t pos, coordf_t gridSize)
{
  float t = (pos / (gridSize * 2.)) + 0.25; // convert relative to grid size
  t = t - (int)t; // extract fractional part
  return((1. - abs(t * 8. - 4.)) * (gridSize / 4.) + (gridSize / 4.));
}

// truncated octagonal waveform, with period and offset
// as per the triangular wave function. The Z position adjusts
// the maximum offset [between -(gridSize / 4) and (gridSize / 4)], with a
// period of (gridSize * 2) and troctWave(Zpos = 0) = 0
static coordf_t troctWave(coordf_t pos, coordf_t gridSize, coordf_t Zpos)
{
  coordf_t Zcycle = triWave(Zpos, gridSize);
  coordf_t perpOffset = Zcycle / 2;
  coordf_t y = triWave(pos, gridSize);
  return((abs(y) > abs(perpOffset)) ?
	 (sgn(y) * perpOffset) :
	 (y * sgn(perpOffset)));
}

// Identify the important points of curve change within a truncated
// octahedron wave (as waveform fraction t):
// 1. Start of wave (always 0.0)
// 2. Transition to upper "horizontal" part
// 3. Transition from upper "horizontal" part
// 4. Transition to lower "horizontal" part
// 5. Transition from lower "horizontal" part
/*    o---o
 *   /     \
 * o/       \
 *           \       /
 *            \     /
 *             o---o
 */
static std::vector<coordf_t> getCriticalPoints(coordf_t Zpos, coordf_t gridSize)
{
  std::vector<coordf_t> res = {0.};
  coordf_t perpOffset = abs(triWave(Zpos, gridSize) / 2.);

  coordf_t normalisedOffset = perpOffset / gridSize;
  // // for debugging: just generate evenly-distributed points
  // for(coordf_t i = 0; i < 2; i += 0.05){
  //   res.push_back(gridSize * i);
  // }
  // note: 0 == straight line
  if(normalisedOffset > 0){
    res.push_back(gridSize * (0. + normalisedOffset));
    res.push_back(gridSize * (1. - normalisedOffset));
    res.push_back(gridSize * (1. + normalisedOffset));
    res.push_back(gridSize * (2. - normalisedOffset));
  }
  return(res);
}

// Generate an array of points that are in the same direction as the
// basic printing line (i.e. Y points for columns, X points for rows)
// Note: a negative offset only causes a change in the perpendicular
// direction
static std::vector<coordf_t> colinearPoints(const coordf_t Zpos, coordf_t gridSize, std::vector<coordf_t> critPoints,
					     const size_t baseLocation, size_t gridLength)
{
  std::vector<coordf_t> points;
  points.push_back(baseLocation);
  for (coordf_t cLoc = baseLocation; cLoc < gridLength; cLoc+= (gridSize*2)) {
    for(size_t pi = 0; pi < critPoints.size(); pi++){
      points.push_back(baseLocation + cLoc + critPoints[pi]);
    }
  }
  points.push_back(gridLength);
  return points;
}

// Generate an array of points for the dimension that is perpendicular to
// the basic printing line (i.e. X points for columns, Y points for rows)
  static std::vector<coordf_t> perpendPoints(const coordf_t Zpos, coordf_t gridSize, std::vector<coordf_t> critPoints,
					     size_t baseLocation, size_t gridLength,
                                             size_t offsetBase, coordf_t perpDir)
{
  std::vector<coordf_t> points;
  points.push_back(offsetBase);
  for (coordf_t cLoc = baseLocation; cLoc < gridLength; cLoc+= gridSize*2) {
    for(size_t pi = 0; pi < critPoints.size(); pi++){
      coordf_t offset = troctWave(critPoints[pi], gridSize, Zpos);
      points.push_back(offsetBase + (offset * perpDir));
    }
  }
  points.push_back(offsetBase);
  return points;
}

static inline Pointfs zip(const std::vector<coordf_t> &x, const std::vector<coordf_t> &y)
{
    assert(x.size() == y.size());
    Pointfs out;
    out.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++ i)
        out.push_back(Vec2d(x[i], y[i]));
    return out;
}

// Generate a set of curves (array of array of 2d points) that describe a
// horizontal slice of a truncated regular octahedron.
static std::vector<Pointfs> makeActualGrid(coordf_t Zpos, coordf_t gridSize, size_t boundsX, size_t boundsY)
{
  std::vector<Pointfs> points;
  std::vector<coordf_t> critPoints = getCriticalPoints(Zpos, gridSize);
  coordf_t zCycle = fmod(Zpos + gridSize/2, gridSize * 2.) / (gridSize * 2.);
  bool printVert = zCycle < 0.5;
  if (printVert) {
    int perpDir = -1;
    for (coordf_t x = 0; x <= (boundsX); x+= gridSize, perpDir *= -1) {
      points.push_back(Pointfs());
      Pointfs &newPoints = points.back();
      newPoints = zip(
		      perpendPoints(Zpos, gridSize, critPoints, 0, boundsY, x, perpDir),
		      colinearPoints(Zpos, gridSize, critPoints, 0, boundsY));
      if (perpDir == 1)
	std::reverse(newPoints.begin(), newPoints.end());
    }
  } else {
    int perpDir = 1;
    for (coordf_t y = gridSize; y <= (boundsY); y+= gridSize, perpDir *= -1) {
      points.push_back(Pointfs());
      Pointfs &newPoints = points.back();
      newPoints = zip(
		      colinearPoints(Zpos, gridSize, critPoints, 0, boundsX),
		      perpendPoints(Zpos, gridSize, critPoints, 0, boundsX, y, perpDir));
      if (perpDir == -1)
	std::reverse(newPoints.begin(), newPoints.end());
    }
  }
  return points;
}

// Generate a set of curves (array of array of 2d points) that describe a
// horizontal slice of a truncated regular octahedron with a specified
// grid square size.
// gridWidth and gridHeight define the width and height of the bounding box respectively
static Polylines makeGrid(coordf_t z, coordf_t gridSize, coordf_t boundWidth, coordf_t boundHeight, bool fillEvenly)
{
  std::vector<Pointfs> polylines = makeActualGrid(z, gridSize, boundWidth, boundHeight);
  Polylines result;
  result.reserve(polylines.size());
  for (std::vector<Pointfs>::const_iterator it_polylines = polylines.begin();
       it_polylines != polylines.end(); ++ it_polylines) {
    result.push_back(Polyline());
    Polyline &polyline = result.back();
    for (Pointfs::const_iterator it = it_polylines->begin(); it != it_polylines->end(); ++ it)
      polyline.points.push_back(Point(coord_t((*it)(0)), coord_t((*it)(1))));
  }
  return result;
}

// FillParams has the following useful information:
// density <0 .. 1>  [proportion of space to fill]
// anchor_length     [???]
// anchor_length_max [???]
// dont_connect()    [avoid connect lines]
// dont_adjust       [avoid filling space evenly]
// monotonic         [fill strictly left to right]
// complete          [complete each loop]
  
void Fill3DHoneycomb::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction, 
    ExPolygon                        expolygon,
    Polylines                       &polylines_out)
{
    // no rotation is supported for this infill pattern
    BoundingBox bb = expolygon.contour.bounding_box();

    // Note: with equally-scaled X/Y/Z, the pattern will create a vertically-stretched
    // truncated octahedron; so Z is pre-adjusted first by scaling by sqrt(2)
    coordf_t zScale = sqrt(2);

    // adjustment to account for the additional distance of octagram curves
    // note: this only strictly applies for a rectangular area where the total
    //       Z travel distance is a multiple of the spacing... but it should
    //       be at least better than the prevous estimate which assumed straight
    //       lines
    // = 4 * integrate(func=4*x(sqrt(2) - 1) + 1, from=0, to=0.25)
    // = (sqrt(2) + 1) / 2 [... I think]
    // make a first guess at the preferred grid Size
    coordf_t gridSize = (scale_(this->spacing) * ((zScale + 1.) / 2.) / params.density);

    // This density calculation is incorrect for many values > 25%, possibly
    // due to quantisation error, so this value is used as a first guess, then the
    // Z scale is adjusted to make the layer patterns consistent / symmetric
    // This means that the resultant infill won't be an ideal truncated octahedron,
    // but it should look better than the equivalent quantised version
    
    coordf_t layerHeight = scale_(thickness_layers);
    // ceiling to an integer value of layers per Z
    // (with a little nudge in case it's close to perfect)
    coordf_t layersPerModule = floor((gridSize * 2) / (zScale * layerHeight) + 0.05);
    if(params.density > 0.42){ // exact layer pattern for >42% density
      layersPerModule = 2;
      // re-adjust the grid size for a partial octahedral path
      // (scale of 1.1 guessed based on modeling)
      gridSize = (scale_(this->spacing) * 1.1 / params.density);
      // re-adjust zScale to make layering consistent
      zScale = (gridSize * 2) / (layersPerModule * layerHeight);
    } else {
      if(layersPerModule < 2){
	layersPerModule = 2;
      }
      // re-adjust zScale to make layering consistent
      zScale = (gridSize * 2) / (layersPerModule * layerHeight);
      // re-adjust the grid size to account for the new zScale
      gridSize = (scale_(this->spacing) * ((zScale + 1.) / 2.) / params.density);
      // re-calculate layersPerModule and zScale
      layersPerModule = floor((gridSize * 2) / (zScale * layerHeight) + 0.05);
      if(layersPerModule < 2){
	layersPerModule = 2;
      }
      zScale = (gridSize * 2) / (layersPerModule * layerHeight);
    }

    // align bounding box to a multiple of our honeycomb grid module
    // (a module is 2*$gridSize since one $gridSize half-module is 
    // growing while the other $gridSize half-module is shrinking)
    bb.merge(align_to_grid(bb.min, Point(gridSize*4, gridSize*4)));
    
    // generate pattern
    Polylines polylines =
      makeGrid(
	       scale_(this->z) * zScale,
	       gridSize,
	       bb.size()(0),
	       bb.size()(1),
	       !params.dont_adjust);

    // Recognize when we're drawing a layer where printVert just flipped.
    // This indicates a transition layer where the 'open squares' of the
    // truncated octahedron faces appear. We fill these squares.
    auto calculate_print_vert = [&](coordf_t z_pos, coordf_t grid_size) {
        // This logic MUST be identical to the one in makeActualGrid to correctly predict its behavior.
        coordf_t z_cycle = fmod(z_pos + grid_size/2., grid_size * 2.) / (grid_size * 2.);
        return z_cycle < 0.5;
    };

    // We must use the single layer height from params to correctly identify the previous layer's Z position.
    // A value <= 0 indicates it's not available, so we skip the check.
    if (params.layer_height > 0) {
        coordf_t single_layer_height_scaled = scale_(params.layer_height);
        coordf_t Zpos_curr = scale_(this->z) * zScale;
        coordf_t Zpos_prev = (scale_(this->z) - single_layer_height_scaled) * zScale;
        coordf_t Zpos_prev_prev = (scale_(this->z) - 2 * single_layer_height_scaled) * zScale;

        bool printVert_curr = calculate_print_vert(Zpos_curr, gridSize);
        bool printVert_prev = calculate_print_vert(Zpos_prev, gridSize);
        bool printVert_prev_prev = calculate_print_vert(Zpos_prev_prev, gridSize);

        if (printVert_curr != printVert_prev) {
            // This is a transition layer. Generate an axial serpentine fill for the open squares.
            Polylines serpentine_polylines;
            serpentine_polylines.reserve((size_t(bb.size()(0) / gridSize) + 2) * (size_t(bb.size()(1) / gridSize) + 2));

            // Determine which checkerboard pattern to use based on the direction of the flip.
            const bool use_even_parity_checkerboard = printVert_curr;

            // 1. Calculate perpOffset to find the real size of the square openings.
            const coordf_t perpOffset_curr = std::abs(triWave(Zpos_curr, gridSize) / 2.0);
            const coordf_t perpOffset_prev = std::abs(triWave(Zpos_prev, gridSize) / 2.0);

            // 2. Calculate the side length of the square.
            const coordf_t square_side_curr = gridSize - 2.0 * perpOffset_curr;
            const coordf_t square_side_prev = gridSize - 2.0 * perpOffset_prev;
            const coordf_t square_side_for_lines = std::max(square_side_curr, square_side_prev);

            // 3. Calculate spacing for the serpentine fill.
            const coordf_t required_spacing = scale_(this->spacing);

            // Only proceed if there is enough space to draw at least one fill line.
            if (square_side_curr > required_spacing) {
                const int num_gaps = ceil(square_side_curr / required_spacing);
                if (num_gaps >= 2) { // Need at least 2 gaps to draw one line in between.
                    const int num_lines = num_gaps + 1;
                    const coordf_t actual_spacing = square_side_curr / static_cast<coordf_t>(num_gaps);

                    const int n_max = ceil(bb.size()(0) / gridSize);
                    const int m_max = ceil(bb.size()(1) / gridSize);

                    for (int n = 0; n <= n_max; ++n) {
                        for (int m = 0; m <= m_max; ++m) {
                            bool is_square_location;
                            if (use_even_parity_checkerboard) {
                                is_square_location = ((n + m) % 2 == 0);
                            } else {
                                is_square_location = ((n + m) % 2 != 0);
                            }

                            if (is_square_location) {
                                // This grid cell corresponds to a small square. Fill it with a serpentine pattern.
                                const coordf_t center_x = n * gridSize + gridSize / 2.;
                                const coordf_t center_y = m * gridSize + gridSize / 2.;
                                
                                Points serpentine_points;
                                // Number of lines to draw is num_lines - 2. Each has 2 points.
                                if (num_lines > 2)
                                    serpentine_points.reserve(2 * (num_lines - 2));

                                if (printVert_curr) { // Vertical (axial) fill pattern
                                    const coordf_t half_side_x = square_side_curr / 2.0;
                                    const coordf_t half_side_y = square_side_for_lines / 2.0;
                                    const coord_t min_x = coord_t(center_x - half_side_x);
                                    const coord_t max_x = coord_t(center_x + half_side_x);
                                    const coord_t min_y = coord_t(center_y - half_side_y);
                                    const coord_t max_y = coord_t(center_y + half_side_y);

                                    for (int i = 1; i <= num_lines - 2; ++i) {
                                        const coord_t current_x = coord_t(min_x + i * actual_spacing);
                                        if (serpentine_points.empty()) {
                                            serpentine_points.emplace_back(current_x, min_y);
                                            serpentine_points.emplace_back(current_x, max_y);
                                        } else {
                                            if (serpentine_points.back().y() == max_y) {
                                                serpentine_points.emplace_back(current_x, max_y);
                                                serpentine_points.emplace_back(current_x, min_y);
                                            } else {
                                                serpentine_points.emplace_back(current_x, min_y);
                                                serpentine_points.emplace_back(current_x, max_y);
                                            }
                                        }
                                    }
                                } else { // Horizontal (axial) fill pattern
                                    const coordf_t half_side_x = square_side_for_lines / 2.0;
                                    const coordf_t half_side_y = square_side_curr / 2.0;
                                    const coord_t min_x = coord_t(center_x - half_side_x);
                                    const coord_t max_x = coord_t(center_x + half_side_x);
                                    const coord_t min_y = coord_t(center_y - half_side_y);
                                    const coord_t max_y = coord_t(center_y + half_side_y);

                                    for (int i = 1; i <= num_lines - 2; ++i) {
                                        const coord_t current_y = coord_t(min_y + i * actual_spacing);
                                        if (serpentine_points.empty()) {
                                            serpentine_points.emplace_back(min_x, current_y);
                                            serpentine_points.emplace_back(max_x, current_y);
                                        } else {
                                            if (serpentine_points.back().x() == max_x) {
                                                serpentine_points.emplace_back(max_x, current_y);
                                                serpentine_points.emplace_back(min_x, current_y);
                                            } else {
                                                serpentine_points.emplace_back(min_x, current_y);
                                                serpentine_points.emplace_back(max_x, current_y);
                                            }
                                        }
                                    }
                                }

                                if (serpentine_points.size() > 1) {
                                    serpentine_polylines.emplace_back(std::move(serpentine_points));
                                }
                            }
                        }
                    }
                }
            }
            // Add the new serpentine polylines to the main pattern.
            append(polylines, std::move(serpentine_polylines));
        } else if ((printVert_prev != printVert_prev_prev) && (printVert_curr == printVert_prev)) {
            // This is the layer AFTER a transition. Draw a transverse serpentine pattern.
            Polylines transverse_polylines;
            
            // Calculate the actual boundaries based on the main infill toolpaths for the current layer.
            const coordf_t perpOffset = std::abs(triWave(Zpos_curr, gridSize) / 2.0);
            const coordf_t hole_side = gridSize - 2.0 * perpOffset;

            const coordf_t scaled_spacing = scale_(this->spacing);
            
            if (hole_side > 0) {
                // The number of lines is based on the total available space in the axial direction.
                const int num_lines = floor(hole_side / scaled_spacing) + 1;

                if (num_lines >= 3) {
                    // The length of the lines is shorter, accounting for switchback clearance.
                    // Each switchback protrudes spacing/2 into the main infill.
                    const coordf_t line_length = hole_side - scaled_spacing;
                    if (line_length <= 0) goto end_transverse_fill; // Skip if lines have no length.

                    // The checkerboard pattern depends on the state of the PREVIOUS layer (the transition layer)
                    const bool use_even_parity_checkerboard = printVert_prev;
                    const int n_max = ceil(bb.size()(0) / gridSize);
                    const int m_max = ceil(bb.size()(1) / gridSize);

                    for (int n = 0; n <= n_max; ++n) {
                        for (int m = 0; m <= m_max; ++m) {
                            bool is_square_location;
                            if (use_even_parity_checkerboard) {
                                is_square_location = ((n + m) % 2 == 0);
                            } else {
                                is_square_location = ((n + m) % 2 != 0);
                            }

                            if (is_square_location) {
                                const coordf_t center_x = n * gridSize + gridSize / 2.;
                                const coordf_t center_y = m * gridSize + gridSize / 2.;
                                
                                Points serpentine_points;
                                serpentine_points.reserve(2 * num_lines);
                                
                                // Center the entire pattern by calculating the offset for the first line.
                                const coordf_t total_pattern_width = (num_lines > 1) ? (num_lines - 1) * scaled_spacing : 0;
                                const coordf_t start_offset = (hole_side - total_pattern_width) / 2.0;

                                if (printVert_curr) { // Main lines are vertical, so draw HORIZONTAL (transverse)
                                    const coordf_t half_line_length = line_length / 2.0;
                                    const coord_t min_x = coord_t(center_x - half_line_length);
                                    const coord_t max_x = coord_t(center_x + half_line_length);
                                    const coordf_t min_y_boundary = center_y - hole_side / 2.0;
                                    
                                    for (int i = 0; i < num_lines; ++i) {
                                        const coord_t current_y = coord_t(min_y_boundary + start_offset + i * scaled_spacing);
                                        if (i % 2 == 0) { // left to right
                                            serpentine_points.emplace_back(min_x, current_y);
                                            serpentine_points.emplace_back(max_x, current_y);
                                        } else { // right to left
                                            serpentine_points.emplace_back(max_x, current_y);
                                            serpentine_points.emplace_back(min_x, current_y);
                                        }
                                    }
                                } else { // Main lines are horizontal, so draw VERTICAL (transverse)
                                    const coordf_t half_line_length = line_length / 2.0;
                                    const coord_t min_y = coord_t(center_y - half_line_length);
                                    const coord_t max_y = coord_t(center_y + half_line_length);
                                    const coordf_t min_x_boundary = center_x - hole_side / 2.0;

                                    for (int i = 0; i < num_lines; ++i) {
                                        const coord_t current_x = coord_t(min_x_boundary + start_offset + i * scaled_spacing);
                                        if (i % 2 == 0) { // top to bottom
                                            serpentine_points.emplace_back(current_x, max_y);
                                            serpentine_points.emplace_back(current_x, min_y);
                                        } else { // bottom to top
                                            serpentine_points.emplace_back(current_x, min_y);
                                            serpentine_points.emplace_back(current_x, max_y);
                                        }
                                    }
                                }
                                if (serpentine_points.size() > 1) {
                                    transverse_polylines.emplace_back(std::move(serpentine_points));
                                }
                            }
                        }
                    }
                }
            }
            end_transverse_fill:;
            append(polylines, std::move(transverse_polylines));
        }
    }

    // move pattern in place
    for (Polyline &pl : polylines){
      pl.translate(bb.min);
    }

    // clip pattern to boundaries, chain the clipped polylines
    polylines = intersection_pl(polylines, expolygon);

    // connect lines if needed
    if (params.dont_connect() || polylines.size() <= 1)
        append(polylines_out, chain_polylines(std::move(polylines)));
    else
        this->connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);
}

} // namespace Slic3r
