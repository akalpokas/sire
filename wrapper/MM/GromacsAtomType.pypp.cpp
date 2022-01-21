// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GromacsAtomType.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/conditional.h"

#include "SireCAS/exp.h"

#include "SireCAS/sum.h"

#include "SireCAS/trigfuncs.h"

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "amberparams.h"

#include "gromacsparams.h"

#include "gromacsparams.h"

SireMM::GromacsAtomType __copy__(const SireMM::GromacsAtomType &other){ return SireMM::GromacsAtomType(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_GromacsAtomType_class(){

    { //::SireMM::GromacsAtomType
        typedef bp::class_< SireMM::GromacsAtomType > GromacsAtomType_exposer_t;
        GromacsAtomType_exposer_t GromacsAtomType_exposer = GromacsAtomType_exposer_t( "GromacsAtomType", "This class represents a Gromacs format atom type. This\ncombines particle type information with charge, mass,\nand LJ parameter info\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope GromacsAtomType_scope( GromacsAtomType_exposer );
        bp::enum_< SireMM::GromacsAtomType::PARTICLE_TYPE>("PARTICLE_TYPE")
            .value("UNKNOWN_TYPE", SireMM::GromacsAtomType::UNKNOWN_TYPE)
            .value("ATOM", SireMM::GromacsAtomType::ATOM)
            .value("SHELL", SireMM::GromacsAtomType::SHELL)
            .value("VIRTUAL", SireMM::GromacsAtomType::VIRTUAL)
            .export_values()
            ;
        GromacsAtomType_exposer.def( bp::init< QString, SireUnits::Dimension::MolarMass >(( bp::arg("atom_type"), bp::arg("mass") ), "Construct, specifying only the mass") );
        GromacsAtomType_exposer.def( bp::init< QString, SireUnits::Dimension::MolarMass, SireUnits::Dimension::Charge, SireMM::GromacsAtomType::PARTICLE_TYPE, SireMM::LJParameter const &, bp::optional< SireMol::Element const & > >(( bp::arg("atom_type"), bp::arg("mass"), bp::arg("charge"), bp::arg("particle_type"), bp::arg("ljparam"), bp::arg("element")=SireMol::Element(0) ), "Construct passing in all parameters") );
        GromacsAtomType_exposer.def( bp::init< QString, QString, SireUnits::Dimension::MolarMass, SireUnits::Dimension::Charge, SireMM::GromacsAtomType::PARTICLE_TYPE, SireMM::LJParameter const &, bp::optional< SireMol::Element const & > >(( bp::arg("atom_type"), bp::arg("bond_type"), bp::arg("mass"), bp::arg("charge"), bp::arg("particle_type"), bp::arg("ljparam"), bp::arg("element")=SireMol::Element(0) ), "Construct passing in all parameters") );
        GromacsAtomType_exposer.def( bp::init< SireMM::GromacsAtomType const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::GromacsAtomType::atomType
        
            typedef ::QString ( ::SireMM::GromacsAtomType::*atomType_function_type)(  ) const;
            atomType_function_type atomType_function_value( &::SireMM::GromacsAtomType::atomType );
            
            GromacsAtomType_exposer.def( 
                "atomType"
                , atomType_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::bondType
        
            typedef ::QString ( ::SireMM::GromacsAtomType::*bondType_function_type)(  ) const;
            bondType_function_type bondType_function_value( &::SireMM::GromacsAtomType::bondType );
            
            GromacsAtomType_exposer.def( 
                "bondType"
                , bondType_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::charge
        
            typedef ::SireUnits::Dimension::Charge ( ::SireMM::GromacsAtomType::*charge_function_type)(  ) const;
            charge_function_type charge_function_value( &::SireMM::GromacsAtomType::charge );
            
            GromacsAtomType_exposer.def( 
                "charge"
                , charge_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::element
        
            typedef ::SireMol::Element ( ::SireMM::GromacsAtomType::*element_function_type)(  ) const;
            element_function_type element_function_value( &::SireMM::GromacsAtomType::element );
            
            GromacsAtomType_exposer.def( 
                "element"
                , element_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::hasMassOnly
        
            typedef bool ( ::SireMM::GromacsAtomType::*hasMassOnly_function_type)(  ) const;
            hasMassOnly_function_type hasMassOnly_function_value( &::SireMM::GromacsAtomType::hasMassOnly );
            
            GromacsAtomType_exposer.def( 
                "hasMassOnly"
                , hasMassOnly_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::isNull
        
            typedef bool ( ::SireMM::GromacsAtomType::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::GromacsAtomType::isNull );
            
            GromacsAtomType_exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::ljParameter
        
            typedef ::SireMM::LJParameter ( ::SireMM::GromacsAtomType::*ljParameter_function_type)(  ) const;
            ljParameter_function_type ljParameter_function_value( &::SireMM::GromacsAtomType::ljParameter );
            
            GromacsAtomType_exposer.def( 
                "ljParameter"
                , ljParameter_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::mass
        
            typedef ::SireUnits::Dimension::MolarMass ( ::SireMM::GromacsAtomType::*mass_function_type)(  ) const;
            mass_function_type mass_function_value( &::SireMM::GromacsAtomType::mass );
            
            GromacsAtomType_exposer.def( 
                "mass"
                , mass_function_value
                , "" );
        
        }
        GromacsAtomType_exposer.def( bp::self != bp::self );
        { //::SireMM::GromacsAtomType::operator=
        
            typedef ::SireMM::GromacsAtomType & ( ::SireMM::GromacsAtomType::*assign_function_type)( ::SireMM::GromacsAtomType const & ) ;
            assign_function_type assign_function_value( &::SireMM::GromacsAtomType::operator= );
            
            GromacsAtomType_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GromacsAtomType_exposer.def( bp::self == bp::self );
        { //::SireMM::GromacsAtomType::particleType
        
            typedef ::SireMM::GromacsAtomType::PARTICLE_TYPE ( ::SireMM::GromacsAtomType::*particleType_function_type)(  ) const;
            particleType_function_type particleType_function_value( &::SireMM::GromacsAtomType::particleType );
            
            GromacsAtomType_exposer.def( 
                "particleType"
                , particleType_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::particleTypeLetter
        
            typedef ::QString ( ::SireMM::GromacsAtomType::*particleTypeLetter_function_type)(  ) const;
            particleTypeLetter_function_type particleTypeLetter_function_value( &::SireMM::GromacsAtomType::particleTypeLetter );
            
            GromacsAtomType_exposer.def( 
                "particleTypeLetter"
                , particleTypeLetter_function_value
                , "Return the single letter that represents the particle type" );
        
        }
        { //::SireMM::GromacsAtomType::particleTypeString
        
            typedef ::QString ( ::SireMM::GromacsAtomType::*particleTypeString_function_type)(  ) const;
            particleTypeString_function_type particleTypeString_function_value( &::SireMM::GromacsAtomType::particleTypeString );
            
            GromacsAtomType_exposer.def( 
                "particleTypeString"
                , particleTypeString_function_value
                , "Return a string version of the particle type" );
        
        }
        { //::SireMM::GromacsAtomType::setElement
        
            typedef void ( ::SireMM::GromacsAtomType::*setElement_function_type)( ::SireMol::Element ) ;
            setElement_function_type setElement_function_value( &::SireMM::GromacsAtomType::setElement );
            
            GromacsAtomType_exposer.def( 
                "setElement"
                , setElement_function_value
                , ( bp::arg("elem") )
                , "Set the element to the passed value. This is useful if the atomtype section has\ninvalid mass informtion, as is the case for many topology files generted by acpype.\nThis allows us to update the element of the type using the mass from the atoms\nsection.\n" );
        
        }
        { //::SireMM::GromacsAtomType::toParticleType
        
            typedef ::SireMM::GromacsAtomType::PARTICLE_TYPE ( *toParticleType_function_type )( ::QString const &,bool * );
            toParticleType_function_type toParticleType_function_value( &::SireMM::GromacsAtomType::toParticleType );
            
            GromacsAtomType_exposer.def( 
                "toParticleType"
                , toParticleType_function_value
                , ( bp::arg("word"), bp::arg("ok")=bp::object() )
                , "Convert the passed string to a Gromacs particle type. Use ok to see if this\nworked correctly" );
        
        }
        { //::SireMM::GromacsAtomType::toString
        
            typedef ::QString ( ::SireMM::GromacsAtomType::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::GromacsAtomType::toString );
            
            GromacsAtomType_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this object" );
        
        }
        { //::SireMM::GromacsAtomType::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::GromacsAtomType::typeName );
            
            GromacsAtomType_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::GromacsAtomType::what
        
            typedef char const * ( ::SireMM::GromacsAtomType::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::GromacsAtomType::what );
            
            GromacsAtomType_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        GromacsAtomType_exposer.staticmethod( "toParticleType" );
        GromacsAtomType_exposer.staticmethod( "typeName" );
        GromacsAtomType_exposer.def( "__copy__", &__copy__);
        GromacsAtomType_exposer.def( "__deepcopy__", &__copy__);
        GromacsAtomType_exposer.def( "clone", &__copy__);
        GromacsAtomType_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::GromacsAtomType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsAtomType_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::GromacsAtomType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsAtomType_exposer.def_pickle(sire_pickle_suite< ::SireMM::GromacsAtomType >());
        GromacsAtomType_exposer.def( "__str__", &__str__< ::SireMM::GromacsAtomType > );
        GromacsAtomType_exposer.def( "__repr__", &__str__< ::SireMM::GromacsAtomType > );
    }

}
