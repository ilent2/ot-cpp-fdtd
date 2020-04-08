/** Include all ot-cpp-fdtd files
 *
 * Not sure if this is a good idea to use this file.
 * I can't remember how some of the classes are implemented and this
 * **might** result in allocating extra static variables.
 * Hopefully the compiler will fix that :P
 *
 * Part of ot-cpp-fdtd, Copyright 2020 Isaac Lenton
 * See LICENSE for details about using/distributing this file.
 */

#include "Simulation.hpp"
#include "Utilities/Vec.hpp"
#include "Utilities/OffsetBox.hpp"
#include "Modules/Indices/Simple.hpp"
#include "Modules/Coordinate/Cartesian.hpp"
#include "Modules/Materials/SimplePermittivity.hpp"
#include "Modules/Materials/HomogeneousPermeability.hpp"
#include "Modules/Materials/TimeAverageH.hpp"
#include "Modules/Materials/TimeAverageE.hpp"
#include "Modules/Tfsf/Simple.hpp"
#include "Sources/Source.hpp"
#include "Sources/PlaneWaveBeam.hpp"
#include "Modules/Cpml/Simple.hpp"
#include "Modules/Output/Report.hpp"
#include "Modules/Output/Integrate/Integrate.hpp"
#include "Modules/Output/Integrate/Output.hpp"
#include "Modules/Output/Integrate/Parameters.hpp"
#include "Modules/Output/Integrate/Iterators.hpp"
#include "Modules/Output/Integrate/Offset.hpp"
#include "Modules/Output/Integrate/MaxwellStressTensor.hpp"
#include "Modules/Output/StrideWindow.hpp"
#include "Modules/Output/Progress.hpp"
#include "Common/Materials.hpp"
#include "Common/Timing.hpp"
#include "Utilities/MakeEvenOdd.hpp"

