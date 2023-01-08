// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "FourAtomFunctions.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/symbols.h"

#include "SireMol/atommatcher.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculeinfodata.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "fouratomfunctions.h"

#include "sireglobal.h"

#include "fouratomfunctions.h"

#include "SireMol/moleculedata.h"

SireMM::FourAtomFunctions __copy__(const SireMM::FourAtomFunctions &other){ return SireMM::FourAtomFunctions(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_FourAtomFunctions_class(){

    { //::SireMM::FourAtomFunctions
        typedef bp::class_< SireMM::FourAtomFunctions, bp::bases< SireMM::AtomFunctions, SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > FourAtomFunctions_exposer_t;
        FourAtomFunctions_exposer_t FourAtomFunctions_exposer = FourAtomFunctions_exposer_t( "FourAtomFunctions", "This class holds the set of FourAtomFunction potentials that\nact between the atoms in a molecule\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope FourAtomFunctions_scope( FourAtomFunctions_exposer );
        FourAtomFunctions_exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") ), "Construct the container to hold the set of four-atom functions\nfor the molecule whose data is in moldata") );
        FourAtomFunctions_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "Construct the container to hold the set of four-atom functions\nfor the molecule whose layout information is in molinfo") );
        FourAtomFunctions_exposer.def( bp::init< SireMM::FourAtomFunctions const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::AtomIdx ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Clear all functions that involve the atom atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::AtomID const & ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Clear any function that acts on the atoms identified by atom\nThrow: SireMol::missing_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )
                , bp::release_gil_policy()
                , "Clear any function that acts between the atoms atom0 to atom2\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )
                , bp::release_gil_policy()
                , "Clear any function that acts between the atoms atom0 to atom3\nThrow: SireMol::missing_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::DihedralID const & ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("dihedralid") )
                , bp::release_gil_policy()
                , "Clear the potential that acts over the dihedral identified by dihedralid\nThis clears all matching dihedrals, so 1-2-3-4 and 4-3-2-1\nThrow: SireMol::missing_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)( ::SireMol::ImproperID const & ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("improperid") )
                , bp::release_gil_policy()
                , "Clear the potential that acts over the improper identified by improperid\nThis clears all matching impropers, so 1-2-3-4 and 1-2-4-3\nThrow: SireMol::missing_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::clear
        
            typedef void ( ::SireMM::FourAtomFunctions::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireMM::FourAtomFunctions::clear );
            
            FourAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , bp::release_gil_policy()
                , "Completely clear all of the functions from this set" );
        
        }
        { //::SireMM::FourAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*force_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::FourAtomFunctions::force );
            
            FourAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the force (derivative of the potential with respect to symbol)\nbetween the atoms atom0 to atom3\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*force_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::FourAtomFunctions::force );
            
            FourAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the force (derivative of the potential with respect to symbol)\nbetween the atoms atom0 to atom3\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*force_function_type)( ::SireMol::DihedralID const &,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::FourAtomFunctions::force );
            
            FourAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("dihedralid"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the force (derivative of the potential with respect to symbol)\non the dihedral identified by dihedralid\nNote that this searches 1-2-3-4 first, then searches for 4-3-2-1\nif no function is found\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*force_function_type)( ::SireMol::ImproperID const &,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::FourAtomFunctions::force );
            
            FourAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("improperid"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the force (derivative of the potential with respect to symbol)\non the improper identified by improperid\nNote that this searches 1-2-3-4 first, then searches for 1-2-4-3\nif no function is found\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::forces
        
            typedef ::QVector< SireMM::FourAtomFunction > ( ::SireMM::FourAtomFunctions::*forces_function_type)( ::SireCAS::Symbol const & ) const;
            forces_function_type forces_function_value( &::SireMM::FourAtomFunctions::forces );
            
            FourAtomFunctions_exposer.def( 
                "forces"
                , forces_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the force functions acting between the identified\nquads of atoms, for the given symbol" );
        
        }
        { //::SireMM::FourAtomFunctions::includeOnly
        
            typedef ::SireMM::FourAtomFunctions ( ::SireMM::FourAtomFunctions::*includeOnly_function_type)( ::SireMol::AtomSelection const &,bool ) const;
            includeOnly_function_type includeOnly_function_value( &::SireMM::FourAtomFunctions::includeOnly );
            
            FourAtomFunctions_exposer.def( 
                "includeOnly"
                , includeOnly_function_value
                , ( bp::arg("selected_atoms"), bp::arg("isstrict")=(bool)(true) )
                , "Return the set of functions where only functions that involve the\natoms in selected_atoms are included. If isstrict is true, then\nonly include functions where all of the atoms are in selected_atoms,\nwhile if isstrict is false, include functions where at least one\natom is in selected_atoms" );
        
        }
        { //::SireMM::FourAtomFunctions::isEmpty
        
            typedef bool ( ::SireMM::FourAtomFunctions::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::FourAtomFunctions::isEmpty );
            
            FourAtomFunctions_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty (has no potentials for any internals)" );
        
        }
        { //::SireMM::FourAtomFunctions::nFunctions
        
            typedef int ( ::SireMM::FourAtomFunctions::*nFunctions_function_type)(  ) const;
            nFunctions_function_type nFunctions_function_value( &::SireMM::FourAtomFunctions::nFunctions );
            
            FourAtomFunctions_exposer.def( 
                "nFunctions"
                , nFunctions_function_value
                , bp::release_gil_policy()
                , "Return the number of functions in this set" );
        
        }
        FourAtomFunctions_exposer.def( bp::self != bp::self );
        { //::SireMM::FourAtomFunctions::operator=
        
            typedef ::SireMM::FourAtomFunctions & ( ::SireMM::FourAtomFunctions::*assign_function_type)( ::SireMM::FourAtomFunctions const & ) ;
            assign_function_type assign_function_value( &::SireMM::FourAtomFunctions::operator= );
            
            FourAtomFunctions_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FourAtomFunctions_exposer.def( bp::self == bp::self );
        { //::SireMM::FourAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*potential_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx ) const;
            potential_function_type potential_function_value( &::SireMM::FourAtomFunctions::potential );
            
            FourAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )
                , bp::release_gil_policy()
                , "Return the function acting between the atoms atom0 to atom3.\nThis returns an empty expression if there is no expression between\nthese atoms\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_atom\n" );
        
        }
        { //::SireMM::FourAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*potential_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            potential_function_type potential_function_value( &::SireMM::FourAtomFunctions::potential );
            
            FourAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") )
                , bp::release_gil_policy()
                , "Return the function acting between the atoms atom0 to atom3.\nThis returns an empty expression if there is no expression between\nthese atoms\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*potential_function_type)( ::SireMol::DihedralID const & ) const;
            potential_function_type potential_function_value( &::SireMM::FourAtomFunctions::potential );
            
            FourAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("dihedralid") )
                , bp::release_gil_policy()
                , "Return the function acting on the dihedral identified by dihedralid.\nThis returns an empty expression if there is no expression on\nthis dihedral\nNote that this searches first for 1-2-3-4, then if no function\nis found, it returns the function for 4-3-2-1\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::FourAtomFunctions::*potential_function_type)( ::SireMol::ImproperID const & ) const;
            potential_function_type potential_function_value( &::SireMM::FourAtomFunctions::potential );
            
            FourAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("improperid") )
                , bp::release_gil_policy()
                , "Return the function acting on the improper identified by improperid.\nThis returns an empty expression if there is no expression on\nthis improper\nNote that this searches first for 1-2-3-4, then if no function\nis found, it returns the function for 1-2-4-3\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::potentials
        
            typedef ::QVector< SireMM::FourAtomFunction > ( ::SireMM::FourAtomFunctions::*potentials_function_type)(  ) const;
            potentials_function_type potentials_function_value( &::SireMM::FourAtomFunctions::potentials );
            
            FourAtomFunctions_exposer.def( 
                "potentials"
                , potentials_function_value
                , bp::release_gil_policy()
                , "Return the potential energy functions acting between the identified\nquads of atoms" );
        
        }
        { //::SireMM::FourAtomFunctions::set
        
            typedef void ( ::SireMM::FourAtomFunctions::*set_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireMol::AtomIdx,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::FourAtomFunctions::set );
            
            FourAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("expression") )
                , bp::release_gil_policy()
                , "Set the potential energy function used by atoms atom0 to atom3\nto be equal to expression - this replaces any existing expression\nIf you want to add an improper then use\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_atom\n" );
        
        }
        { //::SireMM::FourAtomFunctions::set
        
            typedef void ( ::SireMM::FourAtomFunctions::*set_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireMol::AtomID const &,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::FourAtomFunctions::set );
            
            FourAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("expression") )
                , bp::release_gil_policy()
                , "Set the potential energy function used by atoms atom0 to atom3\nto be equal to expression - this replaces any existing expression\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::set
        
            typedef void ( ::SireMM::FourAtomFunctions::*set_function_type)( ::SireMol::DihedralID const &,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::FourAtomFunctions::set );
            
            FourAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("dihedralid"), bp::arg("expression") )
                , bp::release_gil_policy()
                , "Set the potential energy function used for the dihedral identified by dihedralid\nto be equal to expression - this replaces any existing expression\nNote that this will replace any equivalent dihedral, so this will replace\nboth 1-2-3-4 and 4-3-2-1 with the new expression\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::set
        
            typedef void ( ::SireMM::FourAtomFunctions::*set_function_type)( ::SireMol::ImproperID const &,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::FourAtomFunctions::set );
            
            FourAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("improperid"), bp::arg("expression") )
                , bp::release_gil_policy()
                , "Set the potential energy function used for the improper identified by improperid\nto be equal to expression - this replaces any existing expression\nNote that this replaces any existing improper with this function,\ne.g. this replaces both 1-2-3-4 and 1-2-4-3 with this function\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMM::FourAtomFunctions::substitute
        
            typedef void ( ::SireMM::FourAtomFunctions::*substitute_function_type)( ::SireCAS::Identities const & ) ;
            substitute_function_type substitute_function_value( &::SireMM::FourAtomFunctions::substitute );
            
            FourAtomFunctions_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , bp::release_gil_policy()
                , "Perform the substitutions contained in identities in all of\nthe expressions in this set. This could be useful if you have\ndefined these expressions with respect to a lambda parameter,\nand now want to set that value of lambda" );
        
        }
        { //::SireMM::FourAtomFunctions::toString
        
            typedef ::QString ( ::SireMM::FourAtomFunctions::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::FourAtomFunctions::toString );
            
            FourAtomFunctions_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation" );
        
        }
        { //::SireMM::FourAtomFunctions::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::FourAtomFunctions::typeName );
            
            FourAtomFunctions_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FourAtomFunctions_exposer.staticmethod( "typeName" );
        FourAtomFunctions_exposer.def( "__copy__", &__copy__);
        FourAtomFunctions_exposer.def( "__deepcopy__", &__copy__);
        FourAtomFunctions_exposer.def( "clone", &__copy__);
        FourAtomFunctions_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::FourAtomFunctions >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FourAtomFunctions_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::FourAtomFunctions >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FourAtomFunctions_exposer.def_pickle(sire_pickle_suite< ::SireMM::FourAtomFunctions >());
        FourAtomFunctions_exposer.def( "__str__", &__str__< ::SireMM::FourAtomFunctions > );
        FourAtomFunctions_exposer.def( "__repr__", &__str__< ::SireMM::FourAtomFunctions > );
    }

}
