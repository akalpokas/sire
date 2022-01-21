// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "MolEnergyTable.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculeview.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "energytable.h"

#include "energytable.h"

SireFF::MolEnergyTable __copy__(const SireFF::MolEnergyTable &other){ return SireFF::MolEnergyTable(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_MolEnergyTable_class(){

    { //::SireFF::MolEnergyTable
        typedef bp::class_< SireFF::MolEnergyTable > MolEnergyTable_exposer_t;
        MolEnergyTable_exposer_t MolEnergyTable_exposer = MolEnergyTable_exposer_t( "MolEnergyTable", "This class holds the energies of all of the atoms of\nselected CutGroups in a molecule. The MolEnergyTable is used\nto accumulate all of the energies of these atoms during\nan energy evaluation, and also to control which energies are\nevaluated (as only the energies on atoms in selected CutGroups\nare evaluated). This allows you to provide some control over\nthe calculation, e.g. only placing a few protein residues into\nthe energy table, thereby preventing the energy of all atoms\nin a protein from being evaluated if they arent actually\nnecessary.\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope MolEnergyTable_scope( MolEnergyTable_exposer );
        MolEnergyTable_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "Construct a table to hold the forces on all of the CutGroups that\nare selected in molview. The forces are initialised to zero") );
        MolEnergyTable_exposer.def( bp::init< SireFF::MolEnergyTable const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::MolEnergyTable::add
        
            typedef bool ( ::SireFF::MolEnergyTable::*add_function_type)( ::SireMol::CGAtomIdx const &,::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolEnergyTable::add );
            
            MolEnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("cgatomidx"), bp::arg("force") )
                , "Add the force force onto this table. This ignores\nforces calculated for atoms that are in CutGroups that are\nnot in this table - this returns whether or not the\natom is in this table\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireFF::MolEnergyTable::add
        
            typedef bool ( ::SireFF::MolEnergyTable::*add_function_type)( ::SireMol::AtomSelection const &,::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolEnergyTable::add );
            
            MolEnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("selected_atoms"), bp::arg("force") )
                , "Add the force force onto this table for all of the atoms\nin selected_atoms. This ignores forces calculated for atoms\nthat are in CutGroups that are not in this table - this returns whether\nor not any selected atoms are in this table\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolEnergyTable::add
        
            typedef void ( ::SireFF::MolEnergyTable::*add_function_type)( ::SireFF::MolEnergyTable const & ) ;
            add_function_type add_function_value( &::SireFF::MolEnergyTable::add );
            
            MolEnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , "Add the forces contained in other onto this force table. This will only\nadd the forces for CutGroups that are in both tables" );
        
        }
        { //::SireFF::MolEnergyTable::add
        
            typedef void ( ::SireFF::MolEnergyTable::*add_function_type)( ::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::MolEnergyTable::add );
            
            MolEnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("force") )
                , "Add the force force onto all of the atom points in this table" );
        
        }
        { //::SireFF::MolEnergyTable::divide
        
            typedef void ( ::SireFF::MolEnergyTable::*divide_function_type)( double ) ;
            divide_function_type divide_function_value( &::SireFF::MolEnergyTable::divide );
            
            MolEnergyTable_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("value") )
                , "Divide the force at all atom points by value" );
        
        }
        { //::SireFF::MolEnergyTable::initialise
        
            typedef void ( ::SireFF::MolEnergyTable::*initialise_function_type)(  ) ;
            initialise_function_type initialise_function_value( &::SireFF::MolEnergyTable::initialise );
            
            MolEnergyTable_exposer.def( 
                "initialise"
                , initialise_function_value
                , "Initialise this table - this resets all of the forces back to zero" );
        
        }
        { //::SireFF::MolEnergyTable::map
        
            typedef int ( ::SireFF::MolEnergyTable::*map_function_type)( ::SireMol::CGIdx ) const;
            map_function_type map_function_value( &::SireFF::MolEnergyTable::map );
            
            MolEnergyTable_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("cgidx") )
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::molNum
        
            typedef ::SireMol::MolNum ( ::SireFF::MolEnergyTable::*molNum_function_type)(  ) const;
            molNum_function_type molNum_function_value( &::SireFF::MolEnergyTable::molNum );
            
            MolEnergyTable_exposer.def( 
                "molNum"
                , molNum_function_value
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::molUID
        
            typedef ::QUuid const & ( ::SireFF::MolEnergyTable::*molUID_function_type)(  ) const;
            molUID_function_type molUID_function_value( &::SireFF::MolEnergyTable::molUID );
            
            MolEnergyTable_exposer.def( 
                "molUID"
                , molUID_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::multiply
        
            typedef void ( ::SireFF::MolEnergyTable::*multiply_function_type)( double ) ;
            multiply_function_type multiply_function_value( &::SireFF::MolEnergyTable::multiply );
            
            MolEnergyTable_exposer.def( 
                "multiply"
                , multiply_function_value
                , ( bp::arg("value") )
                , "Multiply the force at all atom points by value" );
        
        }
        { //::SireFF::MolEnergyTable::nCutGroups
        
            typedef int ( ::SireFF::MolEnergyTable::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireFF::MolEnergyTable::nCutGroups );
            
            MolEnergyTable_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::nSelectedCutGroups
        
            typedef int ( ::SireFF::MolEnergyTable::*nSelectedCutGroups_function_type)(  ) const;
            nSelectedCutGroups_function_type nSelectedCutGroups_function_value( &::SireFF::MolEnergyTable::nSelectedCutGroups );
            
            MolEnergyTable_exposer.def( 
                "nSelectedCutGroups"
                , nSelectedCutGroups_function_value
                , "" );
        
        }
        MolEnergyTable_exposer.def( bp::self != bp::self );
        MolEnergyTable_exposer.def( bp::self * bp::other< double >() );
        MolEnergyTable_exposer.def( bp::self + bp::self );
        MolEnergyTable_exposer.def( bp::self + bp::other< SireMaths::Vector >() );
        MolEnergyTable_exposer.def( bp::self - bp::self );
        MolEnergyTable_exposer.def( bp::self - bp::other< SireMaths::Vector >() );
        MolEnergyTable_exposer.def( -bp::self );
        MolEnergyTable_exposer.def( bp::self / bp::other< double >() );
        { //::SireFF::MolEnergyTable::operator=
        
            typedef ::SireFF::MolEnergyTable & ( ::SireFF::MolEnergyTable::*assign_function_type)( ::SireFF::MolEnergyTable const & ) ;
            assign_function_type assign_function_value( &::SireFF::MolEnergyTable::operator= );
            
            MolEnergyTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::operator=
        
            typedef ::SireFF::MolEnergyTable & ( ::SireFF::MolEnergyTable::*assign_function_type)( ::SireMaths::Vector const & ) ;
            assign_function_type assign_function_value( &::SireFF::MolEnergyTable::operator= );
            
            MolEnergyTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("force") )
                , bp::return_self< >()
                , "" );
        
        }
        MolEnergyTable_exposer.def( bp::self == bp::self );
        { //::SireFF::MolEnergyTable::selected
        
            typedef bool ( ::SireFF::MolEnergyTable::*selected_function_type)( ::SireMol::CGIdx ) const;
            selected_function_type selected_function_value( &::SireFF::MolEnergyTable::selected );
            
            MolEnergyTable_exposer.def( 
                "selected"
                , selected_function_value
                , ( bp::arg("cgidx") )
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::selectedAll
        
            typedef bool ( ::SireFF::MolEnergyTable::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireFF::MolEnergyTable::selectedAll );
            
            MolEnergyTable_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::setAll
        
            typedef void ( ::SireFF::MolEnergyTable::*setAll_function_type)( ::SireMaths::Vector const & ) ;
            setAll_function_type setAll_function_value( &::SireFF::MolEnergyTable::setAll );
            
            MolEnergyTable_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("force") )
                , "Set all of the forces at the atom points equal to force" );
        
        }
        { //::SireFF::MolEnergyTable::subtract
        
            typedef bool ( ::SireFF::MolEnergyTable::*subtract_function_type)( ::SireMol::CGAtomIdx const &,::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolEnergyTable::subtract );
            
            MolEnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("cgatomidx"), bp::arg("force") )
                , "Subtract the force force from this table. This ignores\nforces calculated for atoms that are in CutGroups that are\nnot in this table - this returns whether or not the\natom is in this table\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireFF::MolEnergyTable::subtract
        
            typedef bool ( ::SireFF::MolEnergyTable::*subtract_function_type)( ::SireMol::AtomSelection const &,::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolEnergyTable::subtract );
            
            MolEnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("selected_atoms"), bp::arg("force") )
                , "Subtract the force force from this table for all of the atoms\nin selected_atoms. This ignores forces calculated for atoms\nthat are in CutGroups that are not in this table - this returns whether\nor not any selected atoms are in this table\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolEnergyTable::subtract
        
            typedef void ( ::SireFF::MolEnergyTable::*subtract_function_type)( ::SireFF::MolEnergyTable const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolEnergyTable::subtract );
            
            MolEnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") )
                , "Subtract the forces contained in other from this force table. This will only\nsubtract the forces for CutGroups that are in both tables" );
        
        }
        { //::SireFF::MolEnergyTable::subtract
        
            typedef void ( ::SireFF::MolEnergyTable::*subtract_function_type)( ::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::MolEnergyTable::subtract );
            
            MolEnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("force") )
                , "Subtract the force force from all of the atom points in this table" );
        
        }
        { //::SireFF::MolEnergyTable::toVector
        
            typedef ::QVector< SireMaths::Vector > ( ::SireFF::MolEnergyTable::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireFF::MolEnergyTable::toVector );
            
            MolEnergyTable_exposer.def( 
                "toVector"
                , toVector_function_value
                , "Return an array of all of the forces on the atoms, in CGAtomIdx order" );
        
        }
        { //::SireFF::MolEnergyTable::toVector
        
            typedef ::QVector< SireMaths::Vector > ( ::SireFF::MolEnergyTable::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireFF::MolEnergyTable::toVector );
            
            MolEnergyTable_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , "Return an array of all of the forces on the atoms selected in selection\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireFF::MolEnergyTable::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::MolEnergyTable::typeName );
            
            MolEnergyTable_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::MolEnergyTable::what
        
            typedef char const * ( ::SireFF::MolEnergyTable::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::MolEnergyTable::what );
            
            MolEnergyTable_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        MolEnergyTable_exposer.staticmethod( "typeName" );
        MolEnergyTable_exposer.def( "__copy__", &__copy__);
        MolEnergyTable_exposer.def( "__deepcopy__", &__copy__);
        MolEnergyTable_exposer.def( "clone", &__copy__);
        MolEnergyTable_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::MolEnergyTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolEnergyTable_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::MolEnergyTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolEnergyTable_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::MolEnergyTable > );
        MolEnergyTable_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::MolEnergyTable > );
        MolEnergyTable_exposer.def( "__str__", &__str__< ::SireFF::MolEnergyTable > );
        MolEnergyTable_exposer.def( "__repr__", &__str__< ::SireFF::MolEnergyTable > );
        MolEnergyTable_exposer.def( "__len__", &__len_size< ::SireFF::MolEnergyTable > );
    }

}
