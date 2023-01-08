// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "BrokenParser.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/booleanproperty.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireFF/ffdetail.h"

#include "SireIO/errors.h"

#include "SireMM/mmdetail.h"

#include "SireMol/core.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireMol/trajectory.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "filetrajectory.h"

#include "moleculeparser.h"

#include "supplementary.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QFile>

#include <QFileInfo>

#include <QMutex>

#include <QTextStream>

#include "moleculeparser.h"

SireIO::BrokenParser __copy__(const SireIO::BrokenParser &other){ return SireIO::BrokenParser(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_BrokenParser_class(){

    { //::SireIO::BrokenParser
        typedef bp::class_< SireIO::BrokenParser, bp::bases< SireIO::MoleculeParser, SireBase::Property > > BrokenParser_exposer_t;
        BrokenParser_exposer_t BrokenParser_exposer = BrokenParser_exposer_t( "BrokenParser", "This is a broken parser, returned when the file cannot be parsed", bp::init< >("") );
        bp::scope BrokenParser_scope( BrokenParser_exposer );
        BrokenParser_exposer.def( bp::init< QString const &, SireBase::PropertyMap const & >(( bp::arg("filename"), bp::arg("map") ), "") );
        BrokenParser_exposer.def( bp::init< QStringList const &, SireBase::PropertyMap const & >(( bp::arg("lines"), bp::arg("map") ), "") );
        BrokenParser_exposer.def( bp::init< SireSystem::System const &, SireBase::PropertyMap const & >(( bp::arg("system"), bp::arg("map") ), "") );
        BrokenParser_exposer.def( bp::init< QString const &, QString const &, QStringList const & >(( bp::arg("filename"), bp::arg("suffix"), bp::arg("errors") ), "") );
        BrokenParser_exposer.def( bp::init< QString const &, QStringList const & >(( bp::arg("filename"), bp::arg("errors") ), "") );
        BrokenParser_exposer.def( bp::init< SireIO::BrokenParser const & >(( bp::arg("other") ), "") );
        { //::SireIO::BrokenParser::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::BrokenParser::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::BrokenParser::construct );
            
            BrokenParser_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed filename" );
        
        }
        { //::SireIO::BrokenParser::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::BrokenParser::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::BrokenParser::construct );
            
            BrokenParser_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed set of lines" );
        
        }
        { //::SireIO::BrokenParser::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::BrokenParser::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::BrokenParser::construct );
            
            BrokenParser_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , bp::release_gil_policy()
                , "Return this parser constructed from the passed SireSystem::System" );
        
        }
        { //::SireIO::BrokenParser::errorReport
        
            typedef ::QString ( ::SireIO::BrokenParser::*errorReport_function_type)(  ) const;
            errorReport_function_type errorReport_function_value( &::SireIO::BrokenParser::errorReport );
            
            BrokenParser_exposer.def( 
                "errorReport"
                , errorReport_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::BrokenParser::formatDescription
        
            typedef ::QString ( ::SireIO::BrokenParser::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::BrokenParser::formatDescription );
            
            BrokenParser_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , bp::release_gil_policy()
                , "Return a description of the file format" );
        
        }
        { //::SireIO::BrokenParser::formatName
        
            typedef ::QString ( ::SireIO::BrokenParser::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::BrokenParser::formatName );
            
            BrokenParser_exposer.def( 
                "formatName"
                , formatName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::BrokenParser::isBroken
        
            typedef bool ( ::SireIO::BrokenParser::*isBroken_function_type)(  ) const;
            isBroken_function_type isBroken_function_value( &::SireIO::BrokenParser::isBroken );
            
            BrokenParser_exposer.def( 
                "isBroken"
                , isBroken_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::BrokenParser::nAtoms
        
            typedef int ( ::SireIO::BrokenParser::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::BrokenParser::nAtoms );
            
            BrokenParser_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BrokenParser_exposer.def( bp::self != bp::self );
        { //::SireIO::BrokenParser::operator=
        
            typedef ::SireIO::BrokenParser & ( ::SireIO::BrokenParser::*assign_function_type)( ::SireIO::BrokenParser const & ) ;
            assign_function_type assign_function_value( &::SireIO::BrokenParser::operator= );
            
            BrokenParser_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("arg0") )
                , bp::return_self< >()
                , "" );
        
        }
        BrokenParser_exposer.def( bp::self == bp::self );
        { //::SireIO::BrokenParser::toSystem
        
            typedef ::SireSystem::System ( ::SireIO::BrokenParser::*toSystem_function_type)( ::SireBase::PropertyMap const & ) const;
            toSystem_function_type toSystem_function_value( &::SireIO::BrokenParser::toSystem );
            
            BrokenParser_exposer.def( 
                "toSystem"
                , toSystem_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireIO::BrokenParser::toSystem
        
            typedef ::SireSystem::System ( ::SireIO::BrokenParser::*toSystem_function_type)( ::SireIO::MoleculeParser const &,::SireBase::PropertyMap const & ) const;
            toSystem_function_type toSystem_function_value( &::SireIO::BrokenParser::toSystem );
            
            BrokenParser_exposer.def( 
                "toSystem"
                , toSystem_function_value
                , ( bp::arg("other"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireIO::BrokenParser::toSystem
        
            typedef ::SireSystem::System ( ::SireIO::BrokenParser::*toSystem_function_type)( ::QList< SireBase::PropPtr< SireIO::MoleculeParser > > const &,::SireBase::PropertyMap const & ) const;
            toSystem_function_type toSystem_function_value( &::SireIO::BrokenParser::toSystem );
            
            BrokenParser_exposer.def( 
                "toSystem"
                , toSystem_function_value
                , ( bp::arg("others"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireIO::BrokenParser::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::BrokenParser::typeName );
            
            BrokenParser_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BrokenParser_exposer.staticmethod( "typeName" );
        BrokenParser_exposer.def( "__copy__", &__copy__);
        BrokenParser_exposer.def( "__deepcopy__", &__copy__);
        BrokenParser_exposer.def( "clone", &__copy__);
        BrokenParser_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::BrokenParser >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BrokenParser_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::BrokenParser >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BrokenParser_exposer.def_pickle(sire_pickle_suite< ::SireIO::BrokenParser >());
        BrokenParser_exposer.def( "__str__", &__str__< ::SireIO::BrokenParser > );
        BrokenParser_exposer.def( "__repr__", &__str__< ::SireIO::BrokenParser > );
    }

}
