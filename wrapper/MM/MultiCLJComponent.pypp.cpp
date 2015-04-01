// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MultiCLJComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "multicljcomponent.h"

#include "multicljcomponent.h"

SireMM::MultiCLJComponent __copy__(const SireMM::MultiCLJComponent &other){ return SireMM::MultiCLJComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_MultiCLJComponent_class(){

    { //::SireMM::MultiCLJComponent
        typedef bp::class_< SireMM::MultiCLJComponent, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > MultiCLJComponent_exposer_t;
        MultiCLJComponent_exposer_t MultiCLJComponent_exposer = MultiCLJComponent_exposer_t( "MultiCLJComponent", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("name")=SireFF::FFName() )) );
        bp::scope MultiCLJComponent_scope( MultiCLJComponent_exposer );
        MultiCLJComponent_exposer.def( bp::init< SireMM::MultiCLJComponent const & >(( bp::arg("other") )) );
        { //::SireMM::MultiCLJComponent::add
        
            typedef int ( ::SireMM::MultiCLJComponent::*add_function_type )( ::QString ) ;
            add_function_type add_function_value( &::SireMM::MultiCLJComponent::add );
            
            MultiCLJComponent_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMM::MultiCLJComponent::changeEnergy
        
            typedef void ( ::SireMM::MultiCLJComponent::*changeEnergy_function_type )( ::SireFF::FF &,::SireMM::MultiCLJEnergy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireMM::MultiCLJComponent::changeEnergy );
            
            MultiCLJComponent_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("cljnrg") ) );
        
        }
        { //::SireMM::MultiCLJComponent::coulomb
        
            typedef ::SireMM::CoulombComponent const & ( ::SireMM::MultiCLJComponent::*coulomb_function_type )(  ) const;
            coulomb_function_type coulomb_function_value( &::SireMM::MultiCLJComponent::coulomb );
            
            MultiCLJComponent_exposer.def( 
                "coulomb"
                , coulomb_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::coulomb
        
            typedef ::SireMM::CoulombComponent const & ( ::SireMM::MultiCLJComponent::*coulomb_function_type )( ::QString ) const;
            coulomb_function_type coulomb_function_value( &::SireMM::MultiCLJComponent::coulomb );
            
            MultiCLJComponent_exposer.def( 
                "coulomb"
                , coulomb_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::count
        
            typedef int ( ::SireMM::MultiCLJComponent::*count_function_type )(  ) const;
            count_function_type count_function_value( &::SireMM::MultiCLJComponent::count );
            
            MultiCLJComponent_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::indexOf
        
            typedef int ( ::SireMM::MultiCLJComponent::*indexOf_function_type )( ::QString ) const;
            indexOf_function_type indexOf_function_value( &::SireMM::MultiCLJComponent::indexOf );
            
            MultiCLJComponent_exposer.def( 
                "indexOf"
                , indexOf_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMM::MultiCLJComponent::keys
        
            typedef ::QStringList ( ::SireMM::MultiCLJComponent::*keys_function_type )(  ) const;
            keys_function_type keys_function_value( &::SireMM::MultiCLJComponent::keys );
            
            MultiCLJComponent_exposer.def( 
                "keys"
                , keys_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::lj
        
            typedef ::SireMM::LJComponent const & ( ::SireMM::MultiCLJComponent::*lj_function_type )(  ) const;
            lj_function_type lj_function_value( &::SireMM::MultiCLJComponent::lj );
            
            MultiCLJComponent_exposer.def( 
                "lj"
                , lj_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::lj
        
            typedef ::SireMM::LJComponent const & ( ::SireMM::MultiCLJComponent::*lj_function_type )( ::QString ) const;
            lj_function_type lj_function_value( &::SireMM::MultiCLJComponent::lj );
            
            MultiCLJComponent_exposer.def( 
                "lj"
                , lj_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::nKeys
        
            typedef int ( ::SireMM::MultiCLJComponent::*nKeys_function_type )(  ) const;
            nKeys_function_type nKeys_function_value( &::SireMM::MultiCLJComponent::nKeys );
            
            MultiCLJComponent_exposer.def( 
                "nKeys"
                , nKeys_function_value );
        
        }
        MultiCLJComponent_exposer.def( bp::self != bp::self );
        { //::SireMM::MultiCLJComponent::operator=
        
            typedef ::SireMM::MultiCLJComponent & ( ::SireMM::MultiCLJComponent::*assign_function_type )( ::SireMM::MultiCLJComponent const & ) ;
            assign_function_type assign_function_value( &::SireMM::MultiCLJComponent::operator= );
            
            MultiCLJComponent_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        MultiCLJComponent_exposer.def( bp::self == bp::self );
        { //::SireMM::MultiCLJComponent::remove
        
            typedef int ( ::SireMM::MultiCLJComponent::*remove_function_type )( ::QString ) ;
            remove_function_type remove_function_value( &::SireMM::MultiCLJComponent::remove );
            
            MultiCLJComponent_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMM::MultiCLJComponent::removeAll
        
            typedef void ( ::SireMM::MultiCLJComponent::*removeAll_function_type )(  ) ;
            removeAll_function_type removeAll_function_value( &::SireMM::MultiCLJComponent::removeAll );
            
            MultiCLJComponent_exposer.def( 
                "removeAll"
                , removeAll_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::rename
        
            typedef ::SireMM::MultiCLJComponent ( ::SireMM::MultiCLJComponent::*rename_function_type )( ::SireFF::FFName const & ) const;
            rename_function_type rename_function_value( &::SireMM::MultiCLJComponent::rename );
            
            MultiCLJComponent_exposer.def( 
                "rename"
                , rename_function_value
                , ( bp::arg("name") ) );
        
        }
        { //::SireMM::MultiCLJComponent::setEnergy
        
            typedef void ( ::SireMM::MultiCLJComponent::*setEnergy_function_type )( ::SireFF::FF &,::SireMM::MultiCLJEnergy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireMM::MultiCLJComponent::setEnergy );
            
            MultiCLJComponent_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("cljnrg") ) );
        
        }
        { //::SireMM::MultiCLJComponent::size
        
            typedef int ( ::SireMM::MultiCLJComponent::*size_function_type )(  ) const;
            size_function_type size_function_value( &::SireMM::MultiCLJComponent::size );
            
            MultiCLJComponent_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::MultiCLJComponent::*symbols_function_type )(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::MultiCLJComponent::symbols );
            
            MultiCLJComponent_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::toString
        
            typedef ::QString ( ::SireMM::MultiCLJComponent::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireMM::MultiCLJComponent::toString );
            
            MultiCLJComponent_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::total
        
            typedef ::SireMM::CLJComponent const & ( ::SireMM::MultiCLJComponent::*total_function_type )(  ) const;
            total_function_type total_function_value( &::SireMM::MultiCLJComponent::total );
            
            MultiCLJComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::total
        
            typedef ::SireMM::CLJComponent const & ( ::SireMM::MultiCLJComponent::*total_function_type )( ::QString ) const;
            total_function_type total_function_value( &::SireMM::MultiCLJComponent::total );
            
            MultiCLJComponent_exposer.def( 
                "total"
                , total_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMM::MultiCLJComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::MultiCLJComponent::typeName );
            
            MultiCLJComponent_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::MultiCLJComponent::what
        
            typedef char const * ( ::SireMM::MultiCLJComponent::*what_function_type )(  ) const;
            what_function_type what_function_value( &::SireMM::MultiCLJComponent::what );
            
            MultiCLJComponent_exposer.def( 
                "what"
                , what_function_value );
        
        }
        MultiCLJComponent_exposer.staticmethod( "typeName" );
        MultiCLJComponent_exposer.def( "__copy__", &__copy__);
        MultiCLJComponent_exposer.def( "__deepcopy__", &__copy__);
        MultiCLJComponent_exposer.def( "clone", &__copy__);
        MultiCLJComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::MultiCLJComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MultiCLJComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::MultiCLJComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MultiCLJComponent_exposer.def( "__str__", &__str__< ::SireMM::MultiCLJComponent > );
        MultiCLJComponent_exposer.def( "__repr__", &__str__< ::SireMM::MultiCLJComponent > );
        MultiCLJComponent_exposer.def( "__len__", &__len_size< ::SireMM::MultiCLJComponent > );
    }

}
