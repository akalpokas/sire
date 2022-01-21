// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "StretchBendComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/ff.h"

#include "SireStream/datastream.h"

#include "internalcomponent.h"

#include "internalcomponent.h"

SireMM::StretchBendComponent __copy__(const SireMM::StretchBendComponent &other){ return SireMM::StretchBendComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_StretchBendComponent_class(){

    { //::SireMM::StretchBendComponent
        typedef bp::class_< SireMM::StretchBendComponent, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > StretchBendComponent_exposer_t;
        StretchBendComponent_exposer_t StretchBendComponent_exposer = StretchBendComponent_exposer_t( "StretchBendComponent", "This class represents a StretchBend component of a forcefield", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("ffname")=SireFF::FFName() ), "Constructor") );
        bp::scope StretchBendComponent_scope( StretchBendComponent_exposer );
        StretchBendComponent_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("symbol") ), "Construct from a symbol\nThrow: SireError::incompatible_error\n") );
        StretchBendComponent_exposer.def( bp::init< SireMM::StretchBendComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::StretchBendComponent::changeEnergy
        
            typedef void ( ::SireMM::StretchBendComponent::*changeEnergy_function_type)( ::SireFF::FF &,::SireMM::StretchBendEnergy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireMM::StretchBendComponent::changeEnergy );
            
            StretchBendComponent_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("sbnrg") )
                , "Change the component of the energy in the forcefield ff\nby delta" );
        
        }
        { //::SireMM::StretchBendComponent::setEnergy
        
            typedef void ( ::SireMM::StretchBendComponent::*setEnergy_function_type)( ::SireFF::FF &,::SireMM::StretchBendEnergy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireMM::StretchBendComponent::setEnergy );
            
            StretchBendComponent_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("sbnrg") )
                , "Set the component of the energy in the forcefield ff\nto be equal to the passed energy" );
        
        }
        { //::SireMM::StretchBendComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::StretchBendComponent::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::StretchBendComponent::symbols );
            
            StretchBendComponent_exposer.def( 
                "symbols"
                , symbols_function_value
                , "" );
        
        }
        { //::SireMM::StretchBendComponent::total
        
            typedef ::SireMM::StretchBendComponent const & ( ::SireMM::StretchBendComponent::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireMM::StretchBendComponent::total );
            
            StretchBendComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::StretchBendComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::StretchBendComponent::typeName );
            
            StretchBendComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::StretchBendComponent::what
        
            typedef char const * ( ::SireMM::StretchBendComponent::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::StretchBendComponent::what );
            
            StretchBendComponent_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        StretchBendComponent_exposer.staticmethod( "typeName" );
        StretchBendComponent_exposer.def( "__copy__", &__copy__);
        StretchBendComponent_exposer.def( "__deepcopy__", &__copy__);
        StretchBendComponent_exposer.def( "clone", &__copy__);
        StretchBendComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::StretchBendComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        StretchBendComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::StretchBendComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        StretchBendComponent_exposer.def( "__setstate__", &__setstate__base64< ::SireMM::StretchBendComponent > );
        StretchBendComponent_exposer.def( "__getstate__", &__getstate__base64< ::SireMM::StretchBendComponent > );
        StretchBendComponent_exposer.def( "__str__", &__str__< ::SireMM::StretchBendComponent > );
        StretchBendComponent_exposer.def( "__repr__", &__str__< ::SireMM::StretchBendComponent > );
        StretchBendComponent_exposer.def( "__hash__", &::SireMM::StretchBendComponent::hash );
    }

}
