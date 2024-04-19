/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_FORCEFIELD_RESPIRATIONCONTROLLER_CPP

#include "respirationController/RespirationController.inl"
#include "sofa/core/ObjectFactory.h"
//#include <sofa/gpu/cuda/CudaFixedConstraint.h>
#include <SofaCUDA/component/constraint/projective/CudaFixedConstraint.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::type;
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(RespirationController)

int CudaRespirationControllerClass = core::RegisterObject("Constant forces applied to given degrees of freedom")
//.add< RespirationController<gpu::cuda::CudaVec3dTypes> >()
.add< RespirationController<gpu::cuda::CudaVec3fTypes> >()
;


} // namespace forcefield

} // namespace component

} // namespace sofa
