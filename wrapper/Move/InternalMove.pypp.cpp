// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "InternalMove.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/angleid.h"

#include "SireMol/atomidx.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/dihedralid.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "flexibility.h"

#include "internalmove.h"

#include <QDebug>

#include <QTime>

#include "internalmove.h"

SireMove::InternalMove __copy__(const SireMove::InternalMove &other){ return SireMove::InternalMove(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_InternalMove_class(){

    { //::SireMove::InternalMove
        typedef bp::class_< SireMove::InternalMove, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > InternalMove_exposer_t;
        InternalMove_exposer_t InternalMove_exposer = InternalMove_exposer_t( "InternalMove", bp::init< bp::optional< SireBase::PropertyMap const & > >(( bp::arg("map")=SireBase::PropertyMap() )) );
        bp::scope InternalMove_scope( InternalMove_exposer );
        InternalMove_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )) );
        InternalMove_exposer.def( bp::init< SireMove::Sampler const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("sampler"), bp::arg("map")=SireBase::PropertyMap() )) );
        InternalMove_exposer.def( bp::init< SireMove::InternalMove const & >(( bp::arg("other") )) );
        { //::SireMove::InternalMove::flexibilityProperty
        
            typedef ::SireBase::PropertyName const & ( ::SireMove::InternalMove::*flexibilityProperty_function_type)(  ) const;
            flexibilityProperty_function_type flexibilityProperty_function_value( &::SireMove::InternalMove::flexibilityProperty );
            
            InternalMove_exposer.def( 
                "flexibilityProperty"
                , flexibilityProperty_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMove::InternalMove::moleculeGroup
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireMove::InternalMove::*moleculeGroup_function_type)(  ) const;
            moleculeGroup_function_type moleculeGroup_function_value( &::SireMove::InternalMove::moleculeGroup );
            
            InternalMove_exposer.def( 
                "moleculeGroup"
                , moleculeGroup_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::InternalMove::move
        
            typedef void ( ::SireMove::InternalMove::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::InternalMove::move );
            
            InternalMove_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) ) );
        
        }
        InternalMove_exposer.def( bp::self != bp::self );
        { //::SireMove::InternalMove::operator=
        
            typedef ::SireMove::InternalMove & ( ::SireMove::InternalMove::*assign_function_type)( ::SireMove::InternalMove const & ) ;
            assign_function_type assign_function_value( &::SireMove::InternalMove::operator= );
            
            InternalMove_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        InternalMove_exposer.def( bp::self == bp::self );
        { //::SireMove::InternalMove::sampler
        
            typedef ::SireMove::Sampler const & ( ::SireMove::InternalMove::*sampler_function_type)(  ) const;
            sampler_function_type sampler_function_value( &::SireMove::InternalMove::sampler );
            
            InternalMove_exposer.def( 
                "sampler"
                , sampler_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireMove::InternalMove::setFlexibilityProperty
        
            typedef void ( ::SireMove::InternalMove::*setFlexibilityProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setFlexibilityProperty_function_type setFlexibilityProperty_function_value( &::SireMove::InternalMove::setFlexibilityProperty );
            
            InternalMove_exposer.def( 
                "setFlexibilityProperty"
                , setFlexibilityProperty_function_value
                , ( bp::arg("property") ) );
        
        }
        { //::SireMove::InternalMove::setGenerator
        
            typedef void ( ::SireMove::InternalMove::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::InternalMove::setGenerator );
            
            InternalMove_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") ) );
        
        }
        { //::SireMove::InternalMove::setSampler
        
            typedef void ( ::SireMove::InternalMove::*setSampler_function_type)( ::SireMove::Sampler const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::InternalMove::setSampler );
            
            InternalMove_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("sampler") ) );
        
        }
        { //::SireMove::InternalMove::setSampler
        
            typedef void ( ::SireMove::InternalMove::*setSampler_function_type)( ::SireMol::MoleculeGroup const & ) ;
            setSampler_function_type setSampler_function_value( &::SireMove::InternalMove::setSampler );
            
            InternalMove_exposer.def( 
                "setSampler"
                , setSampler_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireMove::InternalMove::toString
        
            typedef ::QString ( ::SireMove::InternalMove::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::InternalMove::toString );
            
            InternalMove_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMove::InternalMove::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::InternalMove::typeName );
            
            InternalMove_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        InternalMove_exposer.staticmethod( "typeName" );
        InternalMove_exposer.def( "__copy__", &__copy__);
        InternalMove_exposer.def( "__deepcopy__", &__copy__);
        InternalMove_exposer.def( "clone", &__copy__);
        InternalMove_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::InternalMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalMove_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::InternalMove >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InternalMove_exposer.def( "__str__", &__str__< ::SireMove::InternalMove > );
        InternalMove_exposer.def( "__repr__", &__str__< ::SireMove::InternalMove > );
    }

}
