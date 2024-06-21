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
* This component is not open-source                                           *
*                                                                             *
* Authors: Christian Duriez                                                   *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <Anisotropy/component/initAnisotropy.h>

#include <Anisotropy/component/forcefield/TetrahedronAnisotropicForceField.h>

#include <sofa/simulation/TaskScheduler.h>
#include <sofa/simulation/InitTasks.h>
#include <sofa/helper/AdvancedTimer.h>

#include <SoftRobots.Inverse/component/behavior/Actuator.h>

#include <sofa/simulation/Node.h>

namespace anisotropy::constraint
{

//template <class DataTypes>
//class AnisotropyActuator;

using namespace sofa;
using namespace sofa::defaulttype;
using namespace sofa::type;

using sofa::core::visual::VisualParams ;
using sofa::linearalgebra::BaseVector ;
using sofa::core::ConstraintParams ;
using sofa::core::behavior::Actuator ;
using anisotropy::forcefield::TetrahedronAnisotropicForceField ;

/**
 * This component is used to solve effector constraint by changing young moduli.
 * Description can be found at:
 * https://softrobotscomponents.readthedocs.io
*/
template< class DataTypes >
class AnisotropyControlPointActuator : public Actuator<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AnisotropyControlPointActuator,DataTypes), SOFA_TEMPLATE(Actuator,DataTypes));


//    friend class AnisotropyActuator<DataTypes>;

    typedef typename DataTypes::VecCoord        VecCoord;
    typedef typename DataTypes::VecDeriv        VecDeriv;
    typedef typename DataTypes::VecReal         VecReal;

    typedef typename DataTypes::Coord           Coord;
    typedef typename DataTypes::Deriv           Deriv;
    typedef typename DataTypes::MatrixDeriv     MatrixDeriv;
    typedef typename Coord::value_type          Real;
    typedef Vec<2,Real> Vec2;
    typedef Vec<4,bool> Vec4_bool;
    typedef Vec<5,Real> Vec5;

    typedef typename sofa::core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef Data<VecCoord>		DataVecCoord;
    typedef Data<VecDeriv>		DataVecDeriv;
    typedef Data<MatrixDeriv>    DataMatrixDeriv;


public:
    AnisotropyControlPointActuator(MechanicalState* object = nullptr);
    ~AnisotropyControlPointActuator() override;

    ///////////////// Inherited from BaseObject ///////////////////////////
    void init() override;
    void reinit() override;
    void bwdInit() override;
    void reset() override;
    ///////////////////////////////////////////////////////////////////////


    ///////////////// Inherited from SoftRobotsConstraint /////////////
    void buildConstraintMatrix(const ConstraintParams* cParams ,
                               DataMatrixDeriv &cMatrix,
                               unsigned int &cIndex,
                               const DataVecCoord &x) override;

    void getConstraintViolation(const ConstraintParams* cParams ,
                                BaseVector *resV,
                                const BaseVector *Jdx) override;

    void getConstraintViolation(const ConstraintParams* cParams,
                                BaseVector *resV) override;
    /////////////////////////////////////////////////////////////////////


    ///////////////// Inherited from BaseSoftRobotsConstraint ///////////////////////////
    void storeResults(type::vector<double> &lambda,
                      type::vector<double> &delta) override;
    /////////////////////////////////////////////////////////////////////

protected:

    Data<type::vector<Vec2>> d_boundaries;
    Data<type::vector<Real>> d_maxVariationRatio;
    Data<type::vector<Vec4_bool>> d_anisotropyParameter;
    Data<type::vector<int>> d_cpNumber;
    Data<type::vector<Vec5>> d_targetParamValue;
    Data<SReal >    d_QPobjective;
    Data<bool> d_saveValues;
    Data<std::string> d_nameSlave;

    type::vector<Real> m_initialParamValue;
    type::vector<Real> m_previousParamValue;
    type::vector<Real> m_currentParamValue;
    type::vector<int> m_listCPToWorkOn;
    type::vector<double> m_lambda;

    bool m_initError;
    bool m_isSlave;

//    const DataVecCoord* config_pos;
//    DataVecCoord config_pos_copy;
    vector<VecDeriv> partialDeriv;

    /// Link to be set to the topology container in the component graph.
    TetrahedronAnisotropicForceField< DataTypes > * m_tetraForceField;
    AnisotropyControlPointActuator< DataTypes > * m_masterActuator;
    SingleLink<AnisotropyControlPointActuator<DataTypes>, TetrahedronAnisotropicForceField<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_forceField;
    SingleLink<AnisotropyControlPointActuator<DataTypes>, AnisotropyControlPointActuator<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_masterActuator;
    SingleLink < AnisotropyControlPointActuator<DataTypes>, sofa::simulation::Node , BaseLink::FLAG_STOREPATH > l_nodeToParse;
    MultiLink<AnisotropyControlPointActuator<DataTypes>, AnisotropyControlPointActuator<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_slavesActuator;

    void parseNode(sofa::simulation::Node *node);

    void addParamConstraint(unsigned int constraintId,DataMatrixDeriv &cMatrix,VecDeriv &partialDeriv,bool set=false);

    void saveValuesToFile(const bool append);

private:
    void initLimit();
    void updateLimit();

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using Actuator<DataTypes>::m_state ;
    using Actuator<DataTypes>::m_nbLines ;
    using Actuator<DataTypes>::m_constraintId ;
    using Actuator<DataTypes>::m_hasLambdaMin ;
    using Actuator<DataTypes>::m_hasLambdaMax ;
    using Actuator<DataTypes>::m_lambdaMax ;
    using Actuator<DataTypes>::m_lambdaMin ;
    using Actuator<DataTypes>::getContext ;
    using Actuator<DataTypes>::d_componentState ;
    ///////////////////////////////////////////////////////////////////////////

    class ComputeDfXTask : public simulation::CpuTask
    {
    public:

        ComputeDfXTask(simulation::CpuTask::Status* status): CpuTask(status) {}
        ~ComputeDfXTask() override {}

        MemoryAlloc run() final {
            loopOverCP();
            return MemoryAlloc::Stack;
        }

        void set(AnisotropyControlPointActuator *_actuator)
        {
            actuator = _actuator;
            listCPToWorkOn = actuator->m_listCPToWorkOn;
            anisotropyParameter = actuator->d_anisotropyParameter.getValue();
            tetraForceField = actuator->m_tetraForceField;
//            msg_info("set") << actuator->m_state->read(core::ConstVecCoordId::position())->getValue().size() ;
//            msg_info("set") << actuator->m_nbLines;
            partialDeriv.resize(actuator->m_nbLines);
        }

        void loopOverCP(){
//            sofa::Index i = 0;
            // msg_info("set") << listCPToWorkOn.size() <<","<<anisotropyParameter.size();
//            msg_info("set") << config_pos->getValue().size() <<","<<config_pos;

//            sofa::Index k = 0;
//            for (auto& controlPoint : listCPToWorkOn)
//            {
//                sofa::Index i = 0;
//                for (auto& param : anisotropyParameter[k])
//                {
////                    msg_info("set") << config_pos->getValue().size();
//                    computeDfX(controlPoint,param+1,partialDeriv[i]);
////                    msg_info("loopOverCP") << i;
//                    i++;
//                }
//                k++;
//            }

            sofa::Index k = 0;
            sofa::Index index=0;
            for (auto& controlPoint : listCPToWorkOn)
            {
                for (Index j=0;j<4;j++)
                {
                    if (anisotropyParameter[k][j])
                    {
//                        msg_info("index") << index;
//                        msg_info("index") << j;
                        computeDfX(controlPoint,j+1,partialDeriv[index]);
                        index++;
                    }
                }
                k++;
            }
        }

        void computeDfX(size_t cpNumber,size_t param,VecDeriv &partialDeriv)
        {
            sofa::helper::AdvancedTimer::stepBegin("computeDfX_"+std::to_string(cpNumber)+"_"+std::to_string(param));
            VecDeriv Fxp,Fxpdp;
            type::vector<type::Vec<5,Real>> &controlPoints = *tetraForceField->d_controlPoints.beginEdit();

            Real delta = 0.0001; //0.0001;//controlPoints[cpNumber][param] * 0.000001;

            //////////////////// GET Df(x+d)
//            msg_info("") << "-------------------------- " << delta << "    " << controlPoints[cpNumber][param];
            controlPoints[cpNumber][param] += delta;
            //msg_info() << "-------------------------- " << delta << "    " << controlPoints[cpNumber][param];
            tetraForceField->d_controlPoints.endEdit();
            tetraForceField->reinit();
            getForce(Fxpdp);
            ////////////////////

            //////////////////// GET Df(x)
            controlPoints = *tetraForceField->d_controlPoints.beginEdit();
            controlPoints[cpNumber][param] -= 2 * delta;
//            msg_info("") << "-------------------------- " << delta << "    " << controlPoints[cpNumber][param];
            tetraForceField->d_controlPoints.endEdit();
            tetraForceField->reinit();
            getForce(Fxp);
            ////////////////////

            //////////////////// PUT BACK TO INITIAL VALUE
            controlPoints = *tetraForceField->d_controlPoints.beginEdit();
            controlPoints[cpNumber][param] += delta;
//            msg_info("") << "-------------------------- " << delta << "    " << controlPoints[cpNumber][param];
            tetraForceField->d_controlPoints.endEdit();
            tetraForceField->reinit();
            ////////////////////

            partialDeriv.resize(Fxp.size());
//            msg_info("") << Fxp.size();
            for (unsigned int j=0; j<Fxp.size(); j++)
            {
//                msg_info("") << (Fxpdp[j] - Fxp[j])/(2*delta);
                if (param < 3)
                    partialDeriv[j] = 1000.*(Fxpdp[j] - Fxp[j])/(2*delta);
                else
                    partialDeriv[j] = (Fxpdp[j] - Fxp[j])/(2*delta);
            }
            sofa::helper::AdvancedTimer::stepEnd("computeDfX_"+std::to_string(cpNumber)+"_"+std::to_string(param));
        }

        void getForce(VecDeriv& force)
        {
            DataVecDeriv f;
            VecDeriv f0;

            f0.resize(actuator->m_state->read(core::ConstVecCoordId::position())->getValue().size());

            f.setValue(f0);
            // AddForce(): computes internal forces with respect to given positions and known rest positions.
            //             The velocities v are not used in the computation.
            // void TetrahedronFEMForceField<DataTypes>::addForce (const core::MechanicalParams* /*mparams*/,
            //                                                      DataVecDeriv& d_f,
            //                                                      const DataVecCoord& d_x,
            //                                                      const DataVecDeriv& /*d_v*/)
            DataVecDeriv v;
            tetraForceField->addForce(nullptr, f, actuator->m_state->read(core::ConstVecCoordId::position())->getValue(), v);
            force = f.getValue();
        }

        const DataVecCoord * config_pos = nullptr;
    private:
        AnisotropyControlPointActuator * actuator;
        TetrahedronAnisotropicForceField< DataTypes > * tetraForceField;
        vector<int> listCPToWorkOn;
        vector<Vec4_bool> anisotropyParameter;
        vector<VecDeriv> partialDeriv;
        friend class AnisotropyControlPointActuator;
    };
};


#if !defined(ANISOTROPY_COMPONENT_CONSTRAINT_ANISOTROPYCONTROLPOINTACTUATOR_CPP)
extern template class SOFA_ANISOTROPY_API AnisotropyControlPointActuator<defaulttype::Vec3Types>;
#endif
}
