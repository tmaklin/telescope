// telescope: convert between Themisto and kallisto pseudoalignments
// Copyright (C) 2019 Tommi Mäklin (tommi@maklin.fi)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

#ifndef TELESCOPE_MODE_HPP
#define TELESCOPE_MODE_HPP

#include <string>
#include <exception>
#include <iostream>

namespace telescope {
enum Mode { m_unpaired, m_union, m_intersection };
inline Mode get_mode(const std::string &mode_str) {
    if (mode_str == "unpaired") return m_unpaired;
    if (mode_str == "union") return m_union;
    if (mode_str == "intersection") return m_intersection;
    throw std::runtime_error("Unrecognized paired-end mode.");
}
}

namespace cxxargs {
  inline std::istream& operator>> (std::istream &in, telescope::Mode &t) {
    std::string in_val;
    in >> in_val;
    t = telescope::get_mode(in_val);
    return in;
  }
}

#endif
