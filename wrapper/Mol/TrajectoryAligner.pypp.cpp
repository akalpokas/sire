// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "TrajectoryAligner.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "trajectoryaligner.h"

#include "trajectoryaligner.h"

SireMol::TrajectoryAligner __copy__(const SireMol::TrajectoryAligner &other){ return SireMol::TrajectoryAligner(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_TrajectoryAligner_class(){

    { //::SireMol::TrajectoryAligner
        typedef bp::class_< SireMol::TrajectoryAligner, bp::bases< SireBase::Property > > TrajectoryAligner_exposer_t;
        TrajectoryAligner_exposer_t TrajectoryAligner_exposer = TrajectoryAligner_exposer_t( "TrajectoryAligner", "This class can be used to generate the FrameTransform object\nneeded to align each frame of a trajectory.\n", bp::init< >("") );
        bp::scope TrajectoryAligner_scope( TrajectoryAligner_exposer );
        TrajectoryAligner_exposer.def( bp::init< SireMaths::Vector const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("center"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        TrajectoryAligner_exposer.def( bp::init< SireMol::SelectorM< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        TrajectoryAligner_exposer.def( bp::init< SireMol::TrajectoryAligner const & >(( bp::arg("other") ), "") );
        { //::SireMol::TrajectoryAligner::atoms
        
            typedef ::SireMol::SelectorM< SireMol::Atom > const & ( ::SireMol::TrajectoryAligner::*atoms_function_type)(  ) const;
            atoms_function_type atoms_function_value( &::SireMol::TrajectoryAligner::atoms );
            
            TrajectoryAligner_exposer.def( 
                "atoms"
                , atoms_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::center
        
            typedef ::SireMaths::Vector ( ::SireMol::TrajectoryAligner::*center_function_type)(  ) const;
            center_function_type center_function_value( &::SireMol::TrajectoryAligner::center );
            
            TrajectoryAligner_exposer.def( 
                "center"
                , center_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::count
        
            typedef int ( ::SireMol::TrajectoryAligner::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::TrajectoryAligner::count );
            
            TrajectoryAligner_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::isCompatibleWith
        
            typedef bool ( ::SireMol::TrajectoryAligner::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::TrajectoryAligner::isCompatibleWith );
            
            TrajectoryAligner_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::nFrames
        
            typedef int ( ::SireMol::TrajectoryAligner::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireMol::TrajectoryAligner::nFrames );
            
            TrajectoryAligner_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::nSmooth
        
            typedef int ( ::SireMol::TrajectoryAligner::*nSmooth_function_type)(  ) const;
            nSmooth_function_type nSmooth_function_value( &::SireMol::TrajectoryAligner::nSmooth );
            
            TrajectoryAligner_exposer.def( 
                "nSmooth"
                , nSmooth_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        TrajectoryAligner_exposer.def( bp::self != bp::self );
        { //::SireMol::TrajectoryAligner::operator=
        
            typedef ::SireMol::TrajectoryAligner & ( ::SireMol::TrajectoryAligner::*assign_function_type)( ::SireMol::TrajectoryAligner const & ) ;
            assign_function_type assign_function_value( &::SireMol::TrajectoryAligner::operator= );
            
            TrajectoryAligner_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TrajectoryAligner_exposer.def( bp::self == bp::self );
        { //::SireMol::TrajectoryAligner::operator[]
        
            typedef ::SireMol::FrameTransform ( ::SireMol::TrajectoryAligner::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::TrajectoryAligner::operator[] );
            
            TrajectoryAligner_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::operator[]
        
            typedef ::QList< SireMol::FrameTransform > ( ::SireMol::TrajectoryAligner::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::TrajectoryAligner::operator[] );
            
            TrajectoryAligner_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::operator[]
        
            typedef ::QList< SireMol::FrameTransform > ( ::SireMol::TrajectoryAligner::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::TrajectoryAligner::operator[] );
            
            TrajectoryAligner_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::size
        
            typedef int ( ::SireMol::TrajectoryAligner::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::TrajectoryAligner::size );
            
            TrajectoryAligner_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::toString
        
            typedef ::QString ( ::SireMol::TrajectoryAligner::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::TrajectoryAligner::toString );
            
            TrajectoryAligner_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::TrajectoryAligner::typeName );
            
            TrajectoryAligner_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::what
        
            typedef char const * ( ::SireMol::TrajectoryAligner::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::TrajectoryAligner::what );
            
            TrajectoryAligner_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::TrajectoryAligner::wrap
        
            typedef bool ( ::SireMol::TrajectoryAligner::*wrap_function_type)(  ) const;
            wrap_function_type wrap_function_value( &::SireMol::TrajectoryAligner::wrap );
            
            TrajectoryAligner_exposer.def( 
                "wrap"
                , wrap_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        TrajectoryAligner_exposer.staticmethod( "typeName" );
        TrajectoryAligner_exposer.def( "__copy__", &__copy__);
        TrajectoryAligner_exposer.def( "__deepcopy__", &__copy__);
        TrajectoryAligner_exposer.def( "clone", &__copy__);
        TrajectoryAligner_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::TrajectoryAligner >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrajectoryAligner_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::TrajectoryAligner >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrajectoryAligner_exposer.def_pickle(sire_pickle_suite< ::SireMol::TrajectoryAligner >());
        TrajectoryAligner_exposer.def( "__str__", &__str__< ::SireMol::TrajectoryAligner > );
        TrajectoryAligner_exposer.def( "__repr__", &__str__< ::SireMol::TrajectoryAligner > );
        TrajectoryAligner_exposer.def( "__len__", &__len_size< ::SireMol::TrajectoryAligner > );
    }

}