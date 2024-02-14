// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "LambdaLever.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/values.h"

#include "emle.h"

#include "lambdalever.h"

#include "tostring.h"

#include "SireCAS/values.h"

#include "emle.h"

#include "lambdalever.h"

#include "tostring.h"

SireOpenMM::LambdaLever __copy__(const SireOpenMM::LambdaLever &other){ return SireOpenMM::LambdaLever(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_LambdaLever_class(){

    { //::SireOpenMM::LambdaLever
        typedef bp::class_< SireOpenMM::LambdaLever > LambdaLever_exposer_t;
        LambdaLever_exposer_t LambdaLever_exposer = LambdaLever_exposer_t( "LambdaLever", "This is a lever that is used to change the parameters in an OpenMM\ncontext according to a lambda value. This is actually a collection\nof levers, each of which is controlled by the main lever.\n\nYou can use SireCAS expressions to control how each lever changes\neach parameter\n", bp::init< >("") );
        bp::scope LambdaLever_scope( LambdaLever_exposer );
        LambdaLever_exposer.def( bp::init< SireOpenMM::LambdaLever const & >(( bp::arg("other") ), "") );
        { //::SireOpenMM::LambdaLever::addLever
        
            typedef void ( ::SireOpenMM::LambdaLever::*addLever_function_type)( ::QString const & ) ;
            addLever_function_type addLever_function_value( &::SireOpenMM::LambdaLever::addLever );
            
            LambdaLever_exposer.def( 
                "addLever"
                , addLever_function_value
                , ( bp::arg("lever_name") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::LambdaLever::addPerturbableMolecule
        
            typedef int ( ::SireOpenMM::LambdaLever::*addPerturbableMolecule_function_type)( ::SireOpenMM::OpenMMMolecule const &,::QHash< QString, int > const & ) ;
            addPerturbableMolecule_function_type addPerturbableMolecule_function_value( &::SireOpenMM::LambdaLever::addPerturbableMolecule );
            
            LambdaLever_exposer.def( 
                "addPerturbableMolecule"
                , addPerturbableMolecule_function_value
                , ( bp::arg("molecule"), bp::arg("start_indices") )
                , bp::release_gil_policy()
                , "Add info for the passed perturbable OpenMMMolecule, returning\n  its index in the list of perturbable molecules\n" );
        
        }
        { //::SireOpenMM::LambdaLever::addRestraintIndex
        
            typedef void ( ::SireOpenMM::LambdaLever::*addRestraintIndex_function_type)( ::QString const &,int ) ;
            addRestraintIndex_function_type addRestraintIndex_function_value( &::SireOpenMM::LambdaLever::addRestraintIndex );
            
            LambdaLever_exposer.def( 
                "addRestraintIndex"
                , addRestraintIndex_function_value
                , ( bp::arg("force"), bp::arg("index") )
                , bp::release_gil_policy()
                , "Add the index of a restraint force called restraint in the\n  OpenMM System. There can be multiple restraint forces with\n  the same name\n" );
        
        }
        { //::SireOpenMM::LambdaLever::getForceIndex
        
            typedef int ( ::SireOpenMM::LambdaLever::*getForceIndex_function_type)( ::QString const & ) const;
            getForceIndex_function_type getForceIndex_function_value( &::SireOpenMM::LambdaLever::getForceIndex );
            
            LambdaLever_exposer.def( 
                "getForceIndex"
                , getForceIndex_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Get the index of the force called name. Returns -1 if\n  there is no force with this name\n" );
        
        }
        { //::SireOpenMM::LambdaLever::getForceType
        
            typedef ::QString ( ::SireOpenMM::LambdaLever::*getForceType_function_type)( ::QString const &,::OpenMM::System const & ) const;
            getForceType_function_type getForceType_function_value( &::SireOpenMM::LambdaLever::getForceType );
            
            LambdaLever_exposer.def( 
                "getForceType"
                , getForceType_function_value
                , ( bp::arg("name"), bp::arg("system") )
                , bp::release_gil_policy()
                , "Get the C++ type of the force called name. Returns an\n  empty string if there is no such force\n" );
        
        }
        { //::SireOpenMM::LambdaLever::getPerturbableMoleculeMaps
        
            typedef ::QHash< SireMol::MolNum, SireBase::PropertyMap > ( ::SireOpenMM::LambdaLever::*getPerturbableMoleculeMaps_function_type)(  ) const;
            getPerturbableMoleculeMaps_function_type getPerturbableMoleculeMaps_function_value( &::SireOpenMM::LambdaLever::getPerturbableMoleculeMaps );
            
            LambdaLever_exposer.def( 
                "getPerturbableMoleculeMaps"
                , getPerturbableMoleculeMaps_function_value
                , bp::release_gil_policy()
                , "Return all of the property maps used to find the perturbable properties\n  of the perturbable molecules. This is indexed by molecule number\n" );
        
        }
        { //::SireOpenMM::LambdaLever::getRestraints
        
            typedef ::QList< OpenMM::Force * > ( ::SireOpenMM::LambdaLever::*getRestraints_function_type)( ::QString const &,::OpenMM::System & ) const;
            getRestraints_function_type getRestraints_function_value( &::SireOpenMM::LambdaLever::getRestraints );
            
            LambdaLever_exposer.def( 
                "getRestraints"
                , getRestraints_function_value
                , ( bp::arg("name"), bp::arg("system") )
                , bp::release_gil_policy()
                , "Return the pointers to all of the forces from the passed System\n  are restraints called restraint. This returns an empty\n  list if there are no restraints with this name" );
        
        }
        { //::SireOpenMM::LambdaLever::getSchedule
        
            typedef ::SireCAS::LambdaSchedule ( ::SireOpenMM::LambdaLever::*getSchedule_function_type)(  ) const;
            getSchedule_function_type getSchedule_function_value( &::SireOpenMM::LambdaLever::getSchedule );
            
            LambdaLever_exposer.def( 
                "getSchedule"
                , getSchedule_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::LambdaLever::hasLever
        
            typedef bool ( ::SireOpenMM::LambdaLever::*hasLever_function_type)( ::QString const & ) ;
            hasLever_function_type hasLever_function_value( &::SireOpenMM::LambdaLever::hasLever );
            
            LambdaLever_exposer.def( 
                "hasLever"
                , hasLever_function_value
                , ( bp::arg("lever_name") )
                , bp::release_gil_policy()
                , "" );
        
        }
        LambdaLever_exposer.def( bp::self != bp::self );
        { //::SireOpenMM::LambdaLever::operator=
        
            typedef ::SireOpenMM::LambdaLever & ( ::SireOpenMM::LambdaLever::*assign_function_type)( ::SireOpenMM::LambdaLever const & ) ;
            assign_function_type assign_function_value( &::SireOpenMM::LambdaLever::operator= );
            
            LambdaLever_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        LambdaLever_exposer.def( bp::self == bp::self );
        { //::SireOpenMM::LambdaLever::setConstraintIndicies
        
            typedef void ( ::SireOpenMM::LambdaLever::*setConstraintIndicies_function_type)( int,::QVector< int > const & ) ;
            setConstraintIndicies_function_type setConstraintIndicies_function_value( &::SireOpenMM::LambdaLever::setConstraintIndicies );
            
            LambdaLever_exposer.def( 
                "setConstraintIndicies"
                , setConstraintIndicies_function_value
                , ( bp::arg("idx"), bp::arg("constraint_idxs") )
                , bp::release_gil_policy()
                , "Set the constraint indicies for the perturbable molecule at\n  index mol_idx\n" );
        
        }
        { //::SireOpenMM::LambdaLever::setExceptionIndicies
        
            typedef void ( ::SireOpenMM::LambdaLever::*setExceptionIndicies_function_type)( int,::QString const &,::QVector< boost::tuples::tuple< int, int, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > > const & ) ;
            setExceptionIndicies_function_type setExceptionIndicies_function_value( &::SireOpenMM::LambdaLever::setExceptionIndicies );
            
            LambdaLever_exposer.def( 
                "setExceptionIndicies"
                , setExceptionIndicies_function_value
                , ( bp::arg("idx"), bp::arg("ff"), bp::arg("exception_idxs") )
                , bp::release_gil_policy()
                , "Set the exception indices for the perturbable molecule at\n  index mol_idx\n" );
        
        }
        { //::SireOpenMM::LambdaLever::setForceIndex
        
            typedef void ( ::SireOpenMM::LambdaLever::*setForceIndex_function_type)( ::QString const &,int ) ;
            setForceIndex_function_type setForceIndex_function_value( &::SireOpenMM::LambdaLever::setForceIndex );
            
            LambdaLever_exposer.def( 
                "setForceIndex"
                , setForceIndex_function_value
                , ( bp::arg("force"), bp::arg("index") )
                , bp::release_gil_policy()
                , "Set the index of the force called force in the OpenMM System.\n  There can only be one force with this name. Attempts to add\n  a duplicate will cause an error to be raised.\n" );
        
        }
        { //::SireOpenMM::LambdaLever::setLambda
        
            typedef double ( ::SireOpenMM::LambdaLever::*setLambda_function_type)( ::OpenMM::Context &,double,bool ) const;
            setLambda_function_type setLambda_function_value( &::SireOpenMM::LambdaLever::setLambda );
            
            LambdaLever_exposer.def( 
                "setLambda"
                , setLambda_function_value
                , ( bp::arg("system"), bp::arg("lambda_value"), bp::arg("update_constraints")=(bool)(true) )
                , "Set the value of lambda in the passed context. Returns the\n  actual value of lambda set.\n" );
        
        }
        { //::SireOpenMM::LambdaLever::setSchedule
        
            typedef void ( ::SireOpenMM::LambdaLever::*setSchedule_function_type)( ::SireCAS::LambdaSchedule const & ) ;
            setSchedule_function_type setSchedule_function_value( &::SireOpenMM::LambdaLever::setSchedule );
            
            LambdaLever_exposer.def( 
                "setSchedule"
                , setSchedule_function_value
                , ( bp::arg("schedule") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::LambdaLever::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireOpenMM::LambdaLever::typeName );
            
            LambdaLever_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireOpenMM::LambdaLever::what
        
            typedef char const * ( ::SireOpenMM::LambdaLever::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireOpenMM::LambdaLever::what );
            
            LambdaLever_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        LambdaLever_exposer.staticmethod( "typeName" );
        LambdaLever_exposer.def( "__copy__", &__copy__);
        LambdaLever_exposer.def( "__deepcopy__", &__copy__);
        LambdaLever_exposer.def( "clone", &__copy__);
        LambdaLever_exposer.def( "__str__", &__str__< ::SireOpenMM::LambdaLever > );
        LambdaLever_exposer.def( "__repr__", &__str__< ::SireOpenMM::LambdaLever > );
    }

}
