// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Dynamics.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "dynamics.h"

#include "dynamics.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Dynamics_class(){

    { //::SireMove::Dynamics
        typedef bp::class_< SireMove::Dynamics, bp::bases< SireMove::Move, SireBase::Property >, boost::noncopyable > Dynamics_exposer_t;
        Dynamics_exposer_t Dynamics_exposer = Dynamics_exposer_t( "Dynamics", "This is the base class of the dynamics moves (e.g. all-atom\nmolecular dynamics, or rigid body dynamics)\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Dynamics_scope( Dynamics_exposer );
        { //::SireMove::Dynamics::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::Dynamics::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::Dynamics::ensemble );
            
            Dynamics_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "Return the ensemble that would be generated by this move" );
        
        }
        { //::SireMove::Dynamics::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Dynamics::typeName );
            
            Dynamics_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Dynamics_exposer.staticmethod( "typeName" );
        Dynamics_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Dynamics >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Dynamics_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Dynamics >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Dynamics_exposer.def_pickle(sire_pickle_suite< ::SireMove::Dynamics >());
        Dynamics_exposer.def( "__str__", &__str__< ::SireMove::Dynamics > );
        Dynamics_exposer.def( "__repr__", &__str__< ::SireMove::Dynamics > );
    }

}
