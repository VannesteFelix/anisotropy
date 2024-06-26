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

#include <Anisotropy/component/constraint/AnisotropyControlPointActuator.h>
#include <sofa/helper/cast.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MainTaskSchedulerFactory.h>

#include <iomanip>

namespace anisotropy::constraint
{

using type::vector;
using core::objectmodel::BaseContext;
using sofa::core::objectmodel::ComponentState;

template<class DataTypes>
AnisotropyControlPointActuator<DataTypes>::AnisotropyControlPointActuator(MechanicalState* object)
    : Inherit1(object)

    , d_boundaries(initData(&d_boundaries,"boundaries",""))
    , d_maxVariationRatio(initData(&d_maxVariationRatio, "maxVariationRatio",
                                        "Maximum variation of young / its actual value. \n"
                                        "If unspecified default value 1.0e1."))
    , d_anisotropyParameter(initData(&d_anisotropyParameter,"anisotropyParameter","youngModulusTransverse"))
    , d_cpNumber(initData(&d_cpNumber,"cpNumber","cpNumber"))
    , d_QPobjective(initData(&d_QPobjective,"QPobjective","link to qp objective"))
    , d_targetParamValue(initData(&d_targetParamValue,"targetParamValue","targetParamValue"))
    , d_saveValues(initData(&d_saveValues,false,"saveValues",""))
    , d_nameSlave(initData(&d_nameSlave,"nameSlave",""))
    , m_initError(false)
    , m_isSlave(true)
    , l_forceField(initLink("forceField", "link to the HighOrderTetrahedralCorotationalFEMForceField"))
    , l_masterActuator(initLink("masterActuator","link to actuator master to get constrain index"))
    , l_slavesActuator(initLink("slavesActuator","link to actuator master to get constrain index"))
    , l_nodeToParse(initLink("nodeToParse","link to actuator master to get constrain index"))

{
}


template<class DataTypes>
AnisotropyControlPointActuator<DataTypes>::~AnisotropyControlPointActuator()
{
//    msg_info() << "~AnisotropyControlPointActuator";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::init()
{
//    msg_info() << "ENTER init ";
    d_componentState = ComponentState::Valid;
    Inherit1::init();
    if( not l_masterActuator.empty())
    {
        msg_info() << "I am a slave";
        this->m_constraintType = softrobots::behavior::SoftRobotsBaseConstraint::SLAVE;
    }
    else
    {
        m_isSlave = false;
        if(d_saveValues.getValue())
            saveValuesToFile(false);
    }
//    msg_info() << "EXIT init ";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::reinit()
{
//    msg_info() << "reinit";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::initLimit()
{
//    msg_info() << "ENTER initLimit ";
    m_hasLambdaMin=true;
    m_hasLambdaMax=true;
    updateLimit();
//    msg_info() << "EXIT initLimit ";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::updateLimit()
{
//    msg_info() << "ENTER updateLimit " <<  this->name.getValue() << "   ----------------------------";
//    msg_info() << d_QPobjective.getValue();
//    msg_info() << "lambda min  " <<m_lambdaMin;
//    msg_info() << "lambda max  " <<m_lambdaMax;
//    msg_info() << m_lambdaMin.size();
//    msg_info() << "Current values  " <<m_currentParamValue;
    vector<Real> ratio = d_maxVariationRatio.getValue();
    vector<Vec4_bool> anisotropyParameter = d_anisotropyParameter.getValue();
    vector<Vec2> boundaries = d_boundaries.getValue();
    vector<Vec<5,Real>> controlPoints = m_tetraForceField->d_controlPoints.getValue();

//    msg_info() << d_boundaries.getValue();
//    msg_info() << "d_anisotropyParameter : "<< d_anisotropyParameter.getValue();
    sofa::Index k = 0;
    sofa::Index index=0;
    double norm = 1000.;
    for (auto& controlPoint : m_listCPToWorkOn)
    {
        for (Index j=0;j<4;j++)// auto param : aniParam[k])
        {
            if (anisotropyParameter[k][j])
            {
                m_lambdaMin[index] = 0;
                m_lambdaMax[index] = 0;
                if (j<2)
                {
                    if(anisotropyParameter[k][0] and anisotropyParameter[k][1])
                    {
                        if (j==1)
                        {
        //                    msg_info() << "ENTER updateLimit " <<  this->name.getValue() << "   ----------------------------";
                            /// set Lambda for param 0
                            index -=1;
                            double paramMin=m_currentParamValue[index]/norm-m_currentParamValue[index]/norm*ratio[0];
                            if (paramMin<m_currentParamValue[index+1]/norm)
                            {
        //                        msg_info() << "paramMin : " << paramMin << " " << m_currentParamValue[index+1];
        //                        msg_info() << "m_lambdaMin : "<< m_lambdaMin[index];
//
                                m_lambdaMin[index] = - std::abs(m_currentParamValue[index]/norm - m_currentParamValue[index+1]/norm);
                                if (abs(m_currentParamValue[index] - m_currentParamValue[index+1])>1e-10)
                                    m_lambdaMin[index] = - std::abs(m_currentParamValue[index]/norm - m_currentParamValue[index+1]/norm);
                                else
                                    m_lambdaMin[index] = 0;
                            }
                            else if(paramMin>=boundaries[0][0])
                                m_lambdaMin[index] = -m_currentParamValue[index]/norm*ratio[0];
                            else if(m_currentParamValue[index]/norm - boundaries[0][0]>=0)
                                m_lambdaMin[index] = - std::abs(m_currentParamValue[index]/norm - boundaries[0][0]);

                            double paramMax=m_currentParamValue[index]/norm+m_currentParamValue[index]/norm*ratio[0];
                                if(paramMax<=boundaries[0][1])
                                    m_lambdaMax[index] = m_currentParamValue[index]/norm*ratio[0];
                                else if(boundaries[0][1] - m_currentParamValue[index]/norm >=0)
                                    m_lambdaMax[index] = std::abs(boundaries[0][1] - m_currentParamValue[index]/norm);
                            ////////////////

                            /// set Lambda for param 1
                            index +=1;
                            paramMin=m_currentParamValue[index]/norm-m_currentParamValue[index]/norm*ratio[1];
                //            msg_info() << "---------------------------- : ratio " << ratio[i] << "|  paramMin : "<< paramMin << "  |  "<< -m_currentParamValue[index]*ratio[i];
                            if(paramMin>=boundaries[1][0])
                                //m_lambdaMin=-0.5;
                                m_lambdaMin[index] = -m_currentParamValue[index]/norm*ratio[1];
                            else if(m_currentParamValue[index]/norm - boundaries[1][0]>=0)
                                m_lambdaMin[index] = - std::abs(m_currentParamValue[index]/norm - boundaries[1][0]);
                            paramMax=m_currentParamValue[index]/norm+m_currentParamValue[index]/norm*ratio[1];
                //            msg_info() << "---------------------------- : ratio " << ratio[i] << "|  paramMax : "<< paramMax<< "  |  "<< m_currentParamValue[index]*ratio[i];
                            if (paramMax > m_currentParamValue[index-1]/norm+m_lambdaMin[index-1])
                            {
                                m_lambdaMax[index] = 0;
                                //msg_info() << "paramMax : " << paramMax << " " << m_currentParamValue[index-1]+m_lambdaMin[index-1];
                                //m_lambdaMax[index] = std::abs(m_currentParamValue[index-1] - m_currentParamValue[index]);
                                //msg_info() << "m_lambdaMax : "<< m_lambdaMax[index];
                            }
                            else
                            {
                                if(paramMax<=boundaries[1][1])
                                    m_lambdaMax[index] = m_currentParamValue[index]/norm*ratio[1];
                                else if(boundaries[1][1] - m_currentParamValue[index]/norm >=0)
                                    m_lambdaMax[index] = std::abs(boundaries[1][1] - m_currentParamValue[index]/norm);
                            }
                            ////////////////
                        }
                    }
                    else if (anisotropyParameter[k][0] and not anisotropyParameter[k][1])
                    {
                        double paramMin=m_currentParamValue[index]-m_currentParamValue[index]*ratio[0];

                        if (paramMin<controlPoints[controlPoint][2])
                        {
                            m_lambdaMin[index] = - std::abs(m_currentParamValue[index] - controlPoints[controlPoint][2]);
                        }
                        else if(paramMin>=boundaries[0][0])
                            m_lambdaMin[index] = -m_currentParamValue[index]*ratio[0];
                        else if(m_currentParamValue[index] - boundaries[0][0]>=0)
                            m_lambdaMin[index] = - std::abs(m_currentParamValue[index] - boundaries[0][0]);

                        double paramMax=m_currentParamValue[index]+m_currentParamValue[index]*ratio[0];
                        if(paramMax<=boundaries[0][1])
                            m_lambdaMax[index] = m_currentParamValue[index]*ratio[0];
                        else if(boundaries[0][1] - m_currentParamValue[index] >=0)
                            m_lambdaMax[index] = std::abs(boundaries[0][1] - m_currentParamValue[index]);

                    }
                    else if (not anisotropyParameter[k][0] and anisotropyParameter[k][1])
                    {
                        double paramMin=m_currentParamValue[index]-m_currentParamValue[index]*ratio[1];
                        if(paramMin>=boundaries[1][0])
                            m_lambdaMin[index] = -m_currentParamValue[index]*ratio[1];
                        else if(m_currentParamValue[index] - boundaries[1][0]>=0)
                            m_lambdaMin[index] = - std::abs(m_currentParamValue[index] - boundaries[1][0]);
                        double paramMax=m_currentParamValue[index]+m_currentParamValue[index]*ratio[1];
                        if (paramMax > controlPoints[controlPoint][1])
                        {
                            m_lambdaMax[index] = std::abs(controlPoints[controlPoint][1] - m_currentParamValue[index]);
                        }
                        else
                        {
                            if(paramMax<=boundaries[1][1])
                                m_lambdaMax[index] = m_currentParamValue[index]*ratio[1];
                            else if(boundaries[1][1] - m_currentParamValue[index] >=0)
                                m_lambdaMax[index] = std::abs(boundaries[1][1] - m_currentParamValue[index]);
                        }

                    }
                    else
                    {
                        double paramMin=m_currentParamValue[index]/norm+m_currentParamValue[index]/norm*ratio[j];
            //            msg_info() << "---------------------------- : ratio " << ratio[i] << "|  paramMin : "<< paramMin << "  |  "<< -m_currentParamValue[index]*ratio[i];
                            if(paramMin>=boundaries[j][0])
                                //m_lambdaMin=-0.5;
                                m_lambdaMin[index] = -m_currentParamValue[index]/norm*ratio[j];
//                            else if(m_currentParamValue[index] - boundaries[j][0]>=0)
//                                m_lambdaMin[index] = - std::abs(m_currentParamValue[index] - boundaries[j][0]);
                        double paramMax=m_currentParamValue[index]/norm+m_currentParamValue[index]/norm*ratio[j];
            //            msg_info() << "---------------------------- : ratio " << ratio[i] << "|  paramMax : "<< paramMax<< "  |  "<< m_currentParamValue[index]*ratio[i];
                            if(paramMax<=boundaries[j][1]) // & paramMax<=m_currentParamValue[index-1])
                                //m_lambdaMax= 0.5;
                                m_lambdaMax[index] = m_currentParamValue[index]/norm*ratio[j];
//                            else if(boundaries[j][1] - m_currentParamValue[index] >=0)
//                                m_lambdaMax[index] = std::abs(boundaries[j][1] - m_currentParamValue[index]);
                    }
                }
                else
                {
                    Real delta = ratio[j];
                    double paramMin=m_currentParamValue[index]-delta;
                    if(paramMin>=boundaries[j][0])
                        m_lambdaMin[index] = -delta;
                    else if(m_currentParamValue[index] - boundaries[j][0]>=0)
                        m_lambdaMin[index] = - std::abs(m_currentParamValue[index] - boundaries[j][0]);
                    double paramMax=m_currentParamValue[index]+delta;
                    if(paramMax<=boundaries[j][1])
                        m_lambdaMax[index] =  delta;
                    else if(boundaries[j][1] - m_currentParamValue[index] >=0)
                        m_lambdaMax[index] = std::abs(boundaries[j][1] - m_currentParamValue[index]);
                }

                index++;
            }
        }
        k++;
    }

//    msg_info() << boundaries;
//    msg_info() << "lambda min  " <<m_lambdaMin;
//    msg_info() << "lambda max  " <<m_lambdaMax;
//    msg_info() << "Current values  " <<m_currentParamValue;
//    msg_info() << "EXIT updateLimit " <<  this->name.getValue() << "   ----------------------------";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::bwdInit()
{
    msg_info() << "ENTER bwdInit " <<  this->name.getValue() << "   ----------------------------";

//    if(d_multithreading.getValue())
    // if (not m_isSlave)
    //     simulation::TaskScheduler::getInstance()->init();

    /// Take the user provide topology.
    if (l_forceField.empty())
        msg_error() << "link to Forcefield should be set to ensure right behavior.";

    m_tetraForceField = l_forceField.get();
    msg_info() << "ForceField path used: '" << l_forceField.getLinkedPath() << "'";

    if(m_tetraForceField != nullptr)
    {
        if (not d_cpNumber.getValue().empty()) m_listCPToWorkOn = d_cpNumber.getValue();
        else
        {
            for (size_t j=0;j<m_tetraForceField->d_controlPoints.getValue().size();j++)
                m_listCPToWorkOn.push_back(j);
        }

        if (d_anisotropyParameter.getValue().empty())
        {
            vector<Vec4_bool> anisotropyParameter;

            for (size_t j=0;j<m_listCPToWorkOn.size();j++)
                anisotropyParameter.push_back({true,true,true,true});
            d_anisotropyParameter.setValue(anisotropyParameter);
        }

        const vector<Vec4_bool> &anisotropyParameter = d_anisotropyParameter.getValue();

        m_nbLines=0;
        for (const auto& param : anisotropyParameter)
        {
            for (const auto&i : param)
            {
                if (i)
                    m_nbLines+=1;
            }
        }
//        size_t nb_lines = d_anisotropyParameter.getValue().size()*m_listCPToWorkOn.size();
        msg_info()<< "m_listCPToWorkOn      -----> "<<m_listCPToWorkOn;
        msg_info()<< "d_anisotropyParameter -----> "<<d_anisotropyParameter;
        msg_info()<< "nb_lines              -----> "<<m_nbLines;
        m_lambdaMin.resize(m_nbLines);
        m_lambdaMax.resize(m_nbLines);
        m_initialParamValue.resize(m_nbLines);
        m_currentParamValue.resize(m_nbLines);
        m_previousParamValue.resize(m_nbLines);

        std::fill(m_lambdaMin.begin(), m_lambdaMin.end(), 0);
        std::fill(m_lambdaMax.begin(), m_lambdaMax.end(), 0);
        std::fill(m_initialParamValue.begin(), m_initialParamValue.end(), 0);
        std::fill(m_currentParamValue.begin(), m_currentParamValue.end(), 0);
        std::fill(m_previousParamValue.begin(), m_previousParamValue.end(), 0);

        sofa::Index k = 0;
        sofa::Index index=0;
        for (auto& controlPoint : m_listCPToWorkOn)
        {
            for (Index j=0;j<4;j++)// auto param : aniParam[k])
            {
                if (anisotropyParameter[k][j])
                {
//                    msg_info("index") << index;
//                    msg_info("index") << j;
                    m_currentParamValue[index] = m_tetraForceField->d_controlPoints.getValue()[controlPoint][j+1];
                    m_initialParamValue[index] = m_currentParamValue[index];
                    m_previousParamValue[index] = m_currentParamValue[index];
                    index++;
                }
            }
            k++;
        }
        msg_info()<<m_currentParamValue;
        initLimit();
    }
    else
    {
        msg_error()<<"No TetrahedronFEMForceField found";
        d_componentState = ComponentState::Invalid;
    }


    if(m_isSlave)
    {
        msg_info() << "-----------------------------------> ADD MASTER :" << l_masterActuator.getLinkedPath();
       // m_masterActuator = l_masterActuator.get();
    //        msg_info() << "-----------------------------------> ADD MASTER : "<< m_masterActuator->name;
    }
    else
    if(not m_isSlave)
    {
        simulation::Node* root;
        if( not l_nodeToParse.empty())
            root = l_nodeToParse.get();
        else
            root = down_cast<simulation::Node>(getContext()->getRootContext()->toBaseNode());
        msg_info() << "----------------------------------->  " << root->name.getValue();
        for(auto& child : root->child)
        {
            parseNode(child.get());
        }
    }
    msg_info() << "EXIT bwdInit " <<  this->name.getValue() << "   ----------------------------";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::reset()
{
//    msg_info() << "ENTER reset ";
//    msg_info() << getClass()
//    if(d_componentState.getValue() == ComponentState::Invalid)
//        return;

//    m_currentParamValue = m_initialParamValue;
//    initLimit();
//    msg_info() << "EXIT reset ";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::parseNode(sofa::simulation::Node *node)
{
//    msg_info() << "---------------------> " << node->name.getValue();
//    std::size_t found=node->name.getValue().find(d_nameSlave.getValue());
//    if (found==std::string::npos)
//    {
//        msg_info() << "----  IN  ----";
        for(sofa::core::behavior::BaseConstraintSet * actuator : node->constraintSet)
        {
            if (AnisotropyControlPointActuator *ptr = dynamic_cast<AnisotropyControlPointActuator*>(actuator)){
    //            msg_info() << "-----------------------------------> ACTUATOR FOUND";
                bool toAdd = true;
                for(auto& slave : l_slavesActuator)
                    if (slave->getPathName() == ptr->getPathName())
                        toAdd = false;
                if (toAdd)
                {
                    l_slavesActuator.add(dynamic_cast<AnisotropyControlPointActuator<DataTypes>*>(ptr),ptr->getPathName());
                    ptr->l_masterActuator.set(this);
                    msg_info() << "-----------------------------------> ADD SLAVE :" << actuator->getPathName();
                }
            }
        }

    for(auto& child : node->child){
        if (child != getContext())
            parseNode(child.get());
    }
//    }

}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::addParamConstraint(unsigned int constraintId,
                                                       DataMatrixDeriv &cMatrix,
                                                       VecDeriv &partialDeriv)
{
    MatrixDeriv& matrix = *cMatrix.beginEdit();
    MatrixDerivRowIterator rowIterator = matrix.writeLine(constraintId);

    for (unsigned int j=0; j<partialDeriv.size(); j++)
    {
        if(partialDeriv[j].norm() != 0.0)
        {
            rowIterator.addCol(j, partialDeriv[j]);
        }
    }

    cMatrix.endEdit();
}



template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams,
                                                            DataMatrixDeriv &cMatrix,
                                                            unsigned int &cIndex,
                                                            const DataVecCoord &x)
{
    // msg_info() << "ENTER buildConstraintMatrix for " <<  this->name.getValue() << "   ----------------------------";
    sofa::helper::AdvancedTimer::stepBegin("buildConstraintMatrix_"+this->name.getValue());


    SOFA_UNUSED(cParams);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;


    if (not m_isSlave )
    {
        updateLimit();

        // msg_info() << "MASTER";
        m_constraintId = cIndex;
//        m_nbLines = d_anisotropyParameter.getValue().size()*m_listCPToWorkOn.size();
        cIndex += m_nbLines;

        // simulation::TaskScheduler* taskScheduler = simulation::TaskScheduler::getInstance();
        simulation::TaskScheduler* taskScheduler = simulation::MainTaskSchedulerFactory::createInRegistry();

        simulation::CpuTask::Status status;

        type::vector<AnisotropyControlPointActuator::ComputeDfXTask> tasks;
        sofa::Index nbTasks = l_slavesActuator.size()+1;

        tasks.resize(nbTasks, AnisotropyControlPointActuator::ComputeDfXTask(&status));
        // msg_info() << "task created " << nbTasks;

        tasks[0].set(this);
//        tasks[0].run();
        // msg_info() << "TASK FOR "<< this->name.getValue();
        taskScheduler->addTask(&tasks[0]);

        size_t i = 1;
        for (auto& slave : l_slavesActuator)
        {
            slave->m_nbLines = m_nbLines;
            slave->m_constraintId = m_constraintId;
            tasks[i].set(slave);
             taskScheduler->addTask(&tasks[i]);
//            tasks[i].run();
            i++;

            // msg_info() << "TASK FOR "<< slave->name.getValue();
        }
         taskScheduler->workUntilDone(&status);

        partialDeriv = tasks[0].partialDeriv;
        i = 1;
        for (auto& slave : l_slavesActuator)
        {
            slave->partialDeriv = tasks[i].partialDeriv;
            i++;
        }
    }

    const vector<Vec4_bool> &anisotropyParameter = d_anisotropyParameter.getValue();
    sofa::Index k = 0;
    sofa::Index index=0;
    for (auto& controlPoint : m_listCPToWorkOn)
    {
        for (Index j=0;j<4;j++)
        {
            if (anisotropyParameter[k][j])
            {
                // msg_info() << controlPoint << "  " << m_constraintId << "   " << index;
                addParamConstraint(m_constraintId+index,cMatrix,partialDeriv[index]);
                index++;
            }
        }
        k++;
    }

//    msg_info() << "EXIT buildConstraintMatrix " <<  this->name.getValue() << "   ----------------------------";
    sofa::helper::AdvancedTimer::stepEnd("buildConstraintMatrix_"+this->name.getValue());
}

template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                            BaseVector *resV)
{
    if (not m_isSlave)
        Inherit1::getConstraintViolation(cParams,resV);
}

template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                             BaseVector *resV,const BaseVector *Jdx)
{
//    msg_info() << "--------------------------------------------------------------------------- IN getConstraintViolation";
    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    SOFA_UNUSED(cParams);
//    msg_info()<<m_constraintId << "  " << m_nbLines;
    for(size_t i=0;i<m_nbLines;i++)
    {
        resV->set(m_constraintId+i, 0.0);
    }
//    msg_info() << "--------------------------------------------------------------------------- OUT getConstraintViolation";
}


template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::storeResults(vector<double> &lambda, vector<double> &delta)
{
    // msg_info() << "ENTER storeResults " <<  this->name.getValue() << "   ----------------------------";
    SOFA_UNUSED(delta);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    if (m_isSlave)
        m_lambda = l_masterActuator->m_lambda;
    else
        m_lambda = lambda;

    type::vector<Vec<5,Real>> &controlPoints = *m_tetraForceField->d_controlPoints.beginEdit();
    const vector<Vec4_bool> &anisotropyParameter = d_anisotropyParameter.getValue();

    sofa::Index k = 0;
    sofa::Index index=0;
    for (auto& controlPoint : m_listCPToWorkOn)
    {
        for (Index j=0;j<4;j++)
        {
            if (anisotropyParameter[k][j])
            {
                m_previousParamValue[index] = controlPoints[controlPoint][j+1];
                if (j<2)
                    controlPoints[controlPoint][j+1] += 1000.*m_lambda[index];
                else
                    controlPoints[controlPoint][j+1] += m_lambda[index];
                m_currentParamValue[index] = controlPoints[controlPoint][j+1];
                index++;
            }
        }
        k++;
    }

    m_tetraForceField->d_controlPoints.endEdit();
    m_tetraForceField->reinit();

    if (not m_isSlave)
    {
        std::ostringstream stream;
        stream << std::setiosflags(std::ios::right);
        stream <<"\nCurrent Lambda :\n    " << lambda;
        stream <<"\nWith :\n" <<"    lambda min  " <<m_lambdaMin<<"\n"
                             << "    lambda max  " <<m_lambdaMax<<"\n";
        stream <<"And Boundaries :\n    "<<d_boundaries.getValue()<<"\n";
        stream <<"We optimize on parameters : "<< d_anisotropyParameter.getValue();
        stream <<"\n-----------------------------------";
        if (d_targetParamValue.m_isSet)
        {
            type::vector<Vec5> targetParamValue = d_targetParamValue.getValue();
            bool print_dir = false;
            for (auto& controlPoint : m_listCPToWorkOn)
            {
                for (auto& anisotropyParameter : {0,1})
                    if (targetParamValue[controlPoint][anisotropyParameter+1] != controlPoints[controlPoint][anisotropyParameter+1])
                        stream <<"\nCP"<<controlPoint<<"[YM_"<<anisotropyParameter<<"]"<<":    "<<targetParamValue[controlPoint][anisotropyParameter+1] << " | "
                                       << controlPoints[controlPoint][anisotropyParameter+1] << " | Diff : "
                                       << std::abs(targetParamValue[controlPoint][anisotropyParameter+1]-controlPoints[controlPoint][anisotropyParameter+1]);
                for (auto& anisotropyParameter : {2,3})
                    if (targetParamValue[controlPoint][anisotropyParameter+1] != controlPoints[controlPoint][anisotropyParameter+1])
                    {
                        stream <<"\nCP"<<controlPoint<<"[Angle_"<<anisotropyParameter<<"]"<<": "<<targetParamValue[controlPoint][anisotropyParameter+1]*180.0/M_PI << " | "
                                       << controlPoints[controlPoint][anisotropyParameter+1]*180.0/M_PI << " | Diff : "
                                       << std::abs(targetParamValue[controlPoint][anisotropyParameter+1]*180.0/M_PI-controlPoints[controlPoint][anisotropyParameter+1]*180.0/M_PI);

                        print_dir=true;
                    }
            }
            if (print_dir)
            {
                stream <<"\n-----------------------------------";
                stream <<"\nUnit Vector Direction";
                for (auto& controlPoint : m_listCPToWorkOn)
                {
                    if (targetParamValue[controlPoint][3] != controlPoints[controlPoint][3] or targetParamValue[controlPoint][4] != controlPoints[controlPoint][4])
                        stream <<"\nCP"<<controlPoint<<":\n"
                               << Coord(cos(controlPoints[controlPoint][3])*cos(controlPoints[controlPoint][4]),sin(controlPoints[controlPoint][3])*cos(controlPoints[controlPoint][4]),sin(controlPoints[controlPoint][4]))<< "\n"
                               << Coord(cos(targetParamValue[controlPoint][3])*cos(targetParamValue[controlPoint][4]),sin(targetParamValue[controlPoint][3])*cos(targetParamValue[controlPoint][4]),sin(targetParamValue[controlPoint][4]));
                }
            }
        }
        else
        {
            bool print_dir = false;
            sofa::Index k = 0;
            for (auto& controlPoint : m_listCPToWorkOn)
            {
                for (Index j=0;j<4;j++)
                {
                    if (anisotropyParameter[k][j])
                    {
                        if (j == 0 or j==1)
                                stream <<"\nCP"<<controlPoint<<"[YM_"<<j<<"]:    "<< controlPoints[controlPoint][j+1];
                        if (j == 2 or j==3)
                                stream <<"\nCP"<<controlPoint<<"[Angle_"<<j<<"]: "<< controlPoints[controlPoint][j+1]*180.0/M_PI;
                    }
                }
                k++;
            }
        }

        stream << std::resetiosflags(std::ios::right);
        if (this->f_printLog.getValue())
            msg_info() << stream.str();
    }


    if (not m_isSlave)
    {
//        updateLimit(); // --> update just before computing constrain now
        for(auto& slave : l_slavesActuator)
            slave->storeResults(lambda, delta);
        if(d_saveValues.getValue() and getContext()->isActive())
            saveValuesToFile(true);
    }

    Actuator<DataTypes>::storeResults(lambda, delta);



    // msg_info() << "EXIT storeResults " <<  this->name.getValue() << "   ----------------------------";
}

template<class DataTypes>
void AnisotropyControlPointActuator<DataTypes>::saveValuesToFile(const bool append)
{
    // std::ofstream file;
    // char filename[100];
    // sprintf(filename, "QPInverseProblem_Values.py");
    // if(append)
    // {
    //     file.open(filename, std::ios_base::app);
    //     type::vector<Vec<5,Real>> controlPoints = m_tetraForceField->d_controlPoints.getValue();

    //     const vector<Vec4_bool> &anisotropyParameter = d_anisotropyParameter.getValue();

    //     file << "[";
    //     sofa::Index k = 0;
    //     for (auto& controlPoint : m_listCPToWorkOn)
    //     {
    //         file << "[";
    //         for (Index j=0;j<4;j++)
    //         {
    //             file << controlPoints[controlPoint][j+1];
    //             if (j != 3)
    //                 file << ",";
    //         }
    //         file << "]";
    //         if (controlPoint != m_listCPToWorkOn.back())
    //             file << ",";
    //         k++;
    //     }

    //     file << " ]," << std::endl ;

    // }
    // else
    // {
    //     //file.open(filename);
    //     //file << "[";
    //     file.open(filename, std::ios_base::app);
    //     file << std::endl <<"],[" << std::endl;
    // }


    // file.close();
}

}
