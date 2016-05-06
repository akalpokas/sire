// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Intra14Component.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/ff.h"

#include "SireStream/datastream.h"

#include "internalcomponent.h"

#include "internalcomponent.h"

SireMM::Intra14Component __copy__(const SireMM::Intra14Component &other){ return SireMM::Intra14Component(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Intra14Component_class(){

    { //::SireMM::Intra14Component
        typedef bp::class_< SireMM::Intra14Component, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > Intra14Component_exposer_t;
        Intra14Component_exposer_t Intra14Component_exposer = Intra14Component_exposer_t( "Intra14Component", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("name")=SireFF::FFName() )) );
        bp::scope Intra14Component_scope( Intra14Component_exposer );
        Intra14Component_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("symbol") )) );
        Intra14Component_exposer.def( bp::init< SireMM::Intra14Component const & >(( bp::arg("other") )) );
        { //::SireMM::Intra14Component::changeEnergy
        
            typedef void ( ::SireMM::Intra14Component::*changeEnergy_function_type)( ::SireFF::FF &,::SireMM::Intra14Energy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireMM::Intra14Component::changeEnergy );
            
            Intra14Component_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("nrg") ) );
        
        }
        { //::SireMM::Intra14Component::coulomb
        
            typedef ::SireMM::Intra14CoulombComponent const & ( ::SireMM::Intra14Component::*coulomb_function_type)(  ) const;
            coulomb_function_type coulomb_function_value( &::SireMM::Intra14Component::coulomb );
            
            Intra14Component_exposer.def( 
                "coulomb"
                , coulomb_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::Intra14Component::lj
        
            typedef ::SireMM::Intra14LJComponent const & ( ::SireMM::Intra14Component::*lj_function_type)(  ) const;
            lj_function_type lj_function_value( &::SireMM::Intra14Component::lj );
            
            Intra14Component_exposer.def( 
                "lj"
                , lj_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::Intra14Component::setEnergy
        
            typedef void ( ::SireMM::Intra14Component::*setEnergy_function_type)( ::SireFF::FF &,::SireMM::Intra14Energy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireMM::Intra14Component::setEnergy );
            
            Intra14Component_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("nrg") ) );
        
        }
        { //::SireMM::Intra14Component::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::Intra14Component::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::Intra14Component::symbols );
            
            Intra14Component_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireMM::Intra14Component::total
        
            typedef ::SireMM::Intra14Component const & ( ::SireMM::Intra14Component::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireMM::Intra14Component::total );
            
            Intra14Component_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::Intra14Component::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::Intra14Component::typeName );
            
            Intra14Component_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::Intra14Component::what
        
            typedef char const * ( ::SireMM::Intra14Component::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::Intra14Component::what );
            
            Intra14Component_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Intra14Component_exposer.staticmethod( "typeName" );
        Intra14Component_exposer.def( "__copy__", &__copy__);
        Intra14Component_exposer.def( "__deepcopy__", &__copy__);
        Intra14Component_exposer.def( "clone", &__copy__);
        Intra14Component_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::Intra14Component >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Intra14Component_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::Intra14Component >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Intra14Component_exposer.def( "__str__", &__str__< ::SireMM::Intra14Component > );
        Intra14Component_exposer.def( "__repr__", &__str__< ::SireMM::Intra14Component > );
        Intra14Component_exposer.def( "__hash__", &::SireMM::Intra14Component::hash );
    }

}
