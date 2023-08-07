// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "FrameTransform.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atommatch.h"

#include "trajectoryaligner.h"

#include "trajectoryaligner.h"

SireMol::FrameTransform __copy__(const SireMol::FrameTransform &other){ return SireMol::FrameTransform(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_FrameTransform_class(){

    { //::SireMol::FrameTransform
        typedef bp::class_< SireMol::FrameTransform, bp::bases< SireBase::Property > > FrameTransform_exposer_t;
        FrameTransform_exposer_t FrameTransform_exposer = FrameTransform_exposer_t( "FrameTransform", "This class represents a transformation that is needed to be\nperformed for a specific frame of a trajectory\n", bp::init< >("") );
        bp::scope FrameTransform_scope( FrameTransform_exposer );
        FrameTransform_exposer.def( bp::init< SireBase::PropertyMap const & >(( bp::arg("map") ), "") );
        FrameTransform_exposer.def( bp::init< SireMaths::Transform const &, SireMaths::Vector const &, SireVol::Space const &, int, bool >(( bp::arg("transform"), bp::arg("center"), bp::arg("space"), bp::arg("nsmooth"), bp::arg("wrap") ), "") );
        FrameTransform_exposer.def( bp::init< SireMol::FrameTransform const & >(( bp::arg("other") ), "") );
        { //::SireMol::FrameTransform::apply
        
            typedef ::SireVol::SpacePtr ( ::SireMol::FrameTransform::*apply_function_type)( ::SireVol::Space const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::apply
        
            typedef ::SireMaths::Vector ( ::SireMol::FrameTransform::*apply_function_type)( ::SireMaths::Vector const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::apply
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMol::FrameTransform::*apply_function_type)( ::QVector< SireMaths::Vector > const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::apply
        
            typedef ::SireMaths::Vector ( ::SireMol::FrameTransform::*apply_function_type)( ::SireMaths::Vector const &,::SireVol::Space const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("coords"), bp::arg("space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::apply
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMol::FrameTransform::*apply_function_type)( ::QVector< SireMaths::Vector > const &,::SireVol::Space const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("coords"), bp::arg("space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::apply
        
            typedef ::SireMol::Frame ( ::SireMol::FrameTransform::*apply_function_type)( ::SireMol::Frame const & ) const;
            apply_function_type apply_function_value( &::SireMol::FrameTransform::apply );
            
            FrameTransform_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::center
        
            typedef ::SireMaths::Vector const & ( ::SireMol::FrameTransform::*center_function_type)(  ) const;
            center_function_type center_function_value( &::SireMol::FrameTransform::center );
            
            FrameTransform_exposer.def( 
                "center"
                , center_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::FrameTransform::isCompatibleWith
        
            typedef bool ( ::SireMol::FrameTransform::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::FrameTransform::isCompatibleWith );
            
            FrameTransform_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::isNull
        
            typedef bool ( ::SireMol::FrameTransform::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::FrameTransform::isNull );
            
            FrameTransform_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::nSmooth
        
            typedef int ( ::SireMol::FrameTransform::*nSmooth_function_type)(  ) const;
            nSmooth_function_type nSmooth_function_value( &::SireMol::FrameTransform::nSmooth );
            
            FrameTransform_exposer.def( 
                "nSmooth"
                , nSmooth_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FrameTransform_exposer.def( bp::self != bp::self );
        { //::SireMol::FrameTransform::operator=
        
            typedef ::SireMol::FrameTransform & ( ::SireMol::FrameTransform::*assign_function_type)( ::SireMol::FrameTransform const & ) ;
            assign_function_type assign_function_value( &::SireMol::FrameTransform::operator= );
            
            FrameTransform_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FrameTransform_exposer.def( bp::self == bp::self );
        { //::SireMol::FrameTransform::reverse
        
            typedef ::SireVol::SpacePtr ( ::SireMol::FrameTransform::*reverse_function_type)( ::SireVol::Space const & ) const;
            reverse_function_type reverse_function_value( &::SireMol::FrameTransform::reverse );
            
            FrameTransform_exposer.def( 
                "reverse"
                , reverse_function_value
                , ( bp::arg("space") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::reverse
        
            typedef ::SireMaths::Vector ( ::SireMol::FrameTransform::*reverse_function_type)( ::SireMaths::Vector const & ) const;
            reverse_function_type reverse_function_value( &::SireMol::FrameTransform::reverse );
            
            FrameTransform_exposer.def( 
                "reverse"
                , reverse_function_value
                , ( bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::reverse
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMol::FrameTransform::*reverse_function_type)( ::QVector< SireMaths::Vector > const & ) const;
            reverse_function_type reverse_function_value( &::SireMol::FrameTransform::reverse );
            
            FrameTransform_exposer.def( 
                "reverse"
                , reverse_function_value
                , ( bp::arg("coords") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::reverse
        
            typedef ::SireMol::Frame ( ::SireMol::FrameTransform::*reverse_function_type)( ::SireMol::Frame const & ) const;
            reverse_function_type reverse_function_value( &::SireMol::FrameTransform::reverse );
            
            FrameTransform_exposer.def( 
                "reverse"
                , reverse_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::space
        
            typedef ::SireVol::Space const & ( ::SireMol::FrameTransform::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireMol::FrameTransform::space );
            
            FrameTransform_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMol::FrameTransform::toString
        
            typedef ::QString ( ::SireMol::FrameTransform::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::FrameTransform::toString );
            
            FrameTransform_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::transform
        
            typedef ::SireMaths::Transform const & ( ::SireMol::FrameTransform::*transform_function_type)(  ) const;
            transform_function_type transform_function_value( &::SireMol::FrameTransform::transform );
            
            FrameTransform_exposer.def( 
                "transform"
                , transform_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMol::FrameTransform::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::FrameTransform::typeName );
            
            FrameTransform_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::what
        
            typedef char const * ( ::SireMol::FrameTransform::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::FrameTransform::what );
            
            FrameTransform_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::wouldCreateTransform
        
            typedef bool ( *wouldCreateTransform_function_type )( ::SireBase::PropertyMap const & );
            wouldCreateTransform_function_type wouldCreateTransform_function_value( &::SireMol::FrameTransform::wouldCreateTransform );
            
            FrameTransform_exposer.def( 
                "wouldCreateTransform"
                , wouldCreateTransform_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::FrameTransform::wrap
        
            typedef bool ( ::SireMol::FrameTransform::*wrap_function_type)(  ) const;
            wrap_function_type wrap_function_value( &::SireMol::FrameTransform::wrap );
            
            FrameTransform_exposer.def( 
                "wrap"
                , wrap_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FrameTransform_exposer.staticmethod( "typeName" );
        FrameTransform_exposer.staticmethod( "wouldCreateTransform" );
        FrameTransform_exposer.def( "__copy__", &__copy__);
        FrameTransform_exposer.def( "__deepcopy__", &__copy__);
        FrameTransform_exposer.def( "clone", &__copy__);
        FrameTransform_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::FrameTransform >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FrameTransform_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::FrameTransform >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FrameTransform_exposer.def_pickle(sire_pickle_suite< ::SireMol::FrameTransform >());
        FrameTransform_exposer.def( "__str__", &__str__< ::SireMol::FrameTransform > );
        FrameTransform_exposer.def( "__repr__", &__str__< ::SireMol::FrameTransform > );
    }

}
