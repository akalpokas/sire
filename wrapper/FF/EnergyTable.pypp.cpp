// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "EnergyTable.pypp.hpp"

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

SireFF::EnergyTable __copy__(const SireFF::EnergyTable &other){ return SireFF::EnergyTable(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireFF::EnergyTable&){ return "SireFF::EnergyTable";}

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_EnergyTable_class(){

    { //::SireFF::EnergyTable
        typedef bp::class_< SireFF::EnergyTable > EnergyTable_exposer_t;
        EnergyTable_exposer_t EnergyTable_exposer = EnergyTable_exposer_t( "EnergyTable", "A EnergyTable is a workspace within which all of the energies acting\non the atoms of several molecules may be stored. A EnergyTable is\nused as storing the energies requires lots of memory, and continually\ncreating a deleting such large amouts of memory would be inefficient.\nAlso, using a EnergyTable allows for energies to be accumulated directly,\nrather than requiring intermediate storage space for the\nindividual components.\n\nYou create an energy table to hold all of the energies of all of\nthe atoms of all of the molecules in a specified MoleculeGroup.\nThe energies are held in an array that holds the energies for\nthe molecules in the same order as the molecules appear\nin the molecule group. The energytable also comes with\nan index so you can quickly look up the energies for\na specific molecule.\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope EnergyTable_scope( EnergyTable_exposer );
        EnergyTable_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("molgroup") ), "Construct a table that holds all of the forces on all of the atoms\nfor all of the CutGroups viewed in all of the molecules in the passed\nmolecule group") );
        EnergyTable_exposer.def( bp::init< SireFF::EnergyTable const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::EnergyTable::add
        
            typedef void ( ::SireFF::EnergyTable::*add_function_type)( ::SireFF::EnergyTable const & ) ;
            add_function_type add_function_value( &::SireFF::EnergyTable::add );
            
            EnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Add the contents of the table other onto this table. This will only\nadd the forces for the molecules  grids that are in both tables" );
        
        }
        { //::SireFF::EnergyTable::add
        
            typedef void ( ::SireFF::EnergyTable::*add_function_type)( ::SireMaths::Vector const & ) ;
            add_function_type add_function_value( &::SireFF::EnergyTable::add );
            
            EnergyTable_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("force") )
                , bp::release_gil_policy()
                , "Add the force force onto all of the atom  grid points in this table" );
        
        }
        { //::SireFF::EnergyTable::assertContainsTableFor
        
            typedef void ( ::SireFF::EnergyTable::*assertContainsTableFor_function_type)( ::SireMol::MolNum ) const;
            assertContainsTableFor_function_type assertContainsTableFor_function_value( &::SireFF::EnergyTable::assertContainsTableFor );
            
            EnergyTable_exposer.def( 
                "assertContainsTableFor"
                , assertContainsTableFor_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Assert that this forcetable contains a table for the\nforces for the molecule at number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireFF::EnergyTable::constGetTable
        
            typedef ::SireFF::MolEnergyTable const & ( ::SireFF::EnergyTable::*constGetTable_function_type)( ::SireMol::MolNum ) const;
            constGetTable_function_type constGetTable_function_value( &::SireFF::EnergyTable::constGetTable );
            
            EnergyTable_exposer.def( 
                "constGetTable"
                , constGetTable_function_value
                , ( bp::arg("molnum") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::EnergyTable::containsTable
        
            typedef bool ( ::SireFF::EnergyTable::*containsTable_function_type)( ::SireMol::MolNum ) const;
            containsTable_function_type containsTable_function_value( &::SireFF::EnergyTable::containsTable );
            
            EnergyTable_exposer.def( 
                "containsTable"
                , containsTable_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::EnergyTable::count
        
            typedef int ( ::SireFF::EnergyTable::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireFF::EnergyTable::count );
            
            EnergyTable_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::EnergyTable::divide
        
            typedef void ( ::SireFF::EnergyTable::*divide_function_type)( double ) ;
            divide_function_type divide_function_value( &::SireFF::EnergyTable::divide );
            
            EnergyTable_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Divide the force at all atom and grid points by value" );
        
        }
        { //::SireFF::EnergyTable::getTable
        
            typedef ::SireFF::MolEnergyTable const & ( ::SireFF::EnergyTable::*getTable_function_type)( ::SireMol::MolNum ) const;
            getTable_function_type getTable_function_value( &::SireFF::EnergyTable::getTable );
            
            EnergyTable_exposer.def( 
                "getTable"
                , getTable_function_value
                , ( bp::arg("molnum") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::EnergyTable::index
        
            typedef ::QHash< SireMol::MolNum, int > const & ( ::SireFF::EnergyTable::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireFF::EnergyTable::index );
            
            EnergyTable_exposer.def( 
                "index"
                , index_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::EnergyTable::indexOf
        
            typedef int ( ::SireFF::EnergyTable::*indexOf_function_type)( ::SireMol::MolNum ) const;
            indexOf_function_type indexOf_function_value( &::SireFF::EnergyTable::indexOf );
            
            EnergyTable_exposer.def( 
                "indexOf"
                , indexOf_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Return the index of the molecule with number molnum in this table\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireFF::EnergyTable::initialiseTable
        
            typedef void ( ::SireFF::EnergyTable::*initialiseTable_function_type)( ::SireMol::MolNum ) ;
            initialiseTable_function_type initialiseTable_function_value( &::SireFF::EnergyTable::initialiseTable );
            
            EnergyTable_exposer.def( 
                "initialiseTable"
                , initialiseTable_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Initialise all of the forces for the table for the molecule\nwith number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireFF::EnergyTable::initialiseTables
        
            typedef void ( ::SireFF::EnergyTable::*initialiseTables_function_type)(  ) ;
            initialiseTables_function_type initialiseTables_function_value( &::SireFF::EnergyTable::initialiseTables );
            
            EnergyTable_exposer.def( 
                "initialiseTables"
                , initialiseTables_function_value
                , bp::release_gil_policy()
                , "Initialise all of the tables - this resets all of the forces\nback to zero" );
        
        }
        { //::SireFF::EnergyTable::molNums
        
            typedef ::QList< SireMol::MolNum > ( ::SireFF::EnergyTable::*molNums_function_type)(  ) const;
            molNums_function_type molNums_function_value( &::SireFF::EnergyTable::molNums );
            
            EnergyTable_exposer.def( 
                "molNums"
                , molNums_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::EnergyTable::multiply
        
            typedef void ( ::SireFF::EnergyTable::*multiply_function_type)( double ) ;
            multiply_function_type multiply_function_value( &::SireFF::EnergyTable::multiply );
            
            EnergyTable_exposer.def( 
                "multiply"
                , multiply_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Multiply the force at all atom and grid points by value" );
        
        }
        EnergyTable_exposer.def( bp::self != bp::self );
        EnergyTable_exposer.def( bp::self * bp::other< double >() );
        EnergyTable_exposer.def( bp::self + bp::self );
        EnergyTable_exposer.def( bp::self + bp::other< SireMaths::Vector >() );
        EnergyTable_exposer.def( bp::self - bp::self );
        EnergyTable_exposer.def( bp::self - bp::other< SireMaths::Vector >() );
        EnergyTable_exposer.def( -bp::self );
        EnergyTable_exposer.def( bp::self / bp::other< double >() );
        { //::SireFF::EnergyTable::operator=
        
            typedef ::SireFF::EnergyTable & ( ::SireFF::EnergyTable::*assign_function_type)( ::SireFF::EnergyTable const & ) ;
            assign_function_type assign_function_value( &::SireFF::EnergyTable::operator= );
            
            EnergyTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::EnergyTable::operator=
        
            typedef ::SireFF::EnergyTable & ( ::SireFF::EnergyTable::*assign_function_type)( ::SireMaths::Vector const & ) ;
            assign_function_type assign_function_value( &::SireFF::EnergyTable::operator= );
            
            EnergyTable_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("force") )
                , bp::return_self< >()
                , "" );
        
        }
        EnergyTable_exposer.def( bp::self == bp::self );
        { //::SireFF::EnergyTable::setAll
        
            typedef void ( ::SireFF::EnergyTable::*setAll_function_type)( ::SireMaths::Vector const & ) ;
            setAll_function_type setAll_function_value( &::SireFF::EnergyTable::setAll );
            
            EnergyTable_exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("force") )
                , bp::release_gil_policy()
                , "Set the force at all atom and grid points equal to force" );
        
        }
        { //::SireFF::EnergyTable::subtract
        
            typedef void ( ::SireFF::EnergyTable::*subtract_function_type)( ::SireFF::EnergyTable const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::EnergyTable::subtract );
            
            EnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Subtract the contents of the table other from this table. This will only\nsubtract the forces for the molecules  grids that are in both tables" );
        
        }
        { //::SireFF::EnergyTable::subtract
        
            typedef void ( ::SireFF::EnergyTable::*subtract_function_type)( ::SireMaths::Vector const & ) ;
            subtract_function_type subtract_function_value( &::SireFF::EnergyTable::subtract );
            
            EnergyTable_exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("force") )
                , bp::release_gil_policy()
                , "Subtract the force force from all of the atom  grid points in this table" );
        
        }
        { //::SireFF::EnergyTable::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::EnergyTable::typeName );
            
            EnergyTable_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::EnergyTable::what
        
            typedef char const * ( ::SireFF::EnergyTable::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::EnergyTable::what );
            
            EnergyTable_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        EnergyTable_exposer.staticmethod( "typeName" );
        EnergyTable_exposer.def( "__copy__", &__copy__);
        EnergyTable_exposer.def( "__deepcopy__", &__copy__);
        EnergyTable_exposer.def( "clone", &__copy__);
        EnergyTable_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::EnergyTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EnergyTable_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::EnergyTable >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EnergyTable_exposer.def_pickle(sire_pickle_suite< ::SireFF::EnergyTable >());
        EnergyTable_exposer.def( "__str__", &pvt_get_name);
        EnergyTable_exposer.def( "__repr__", &pvt_get_name);
        EnergyTable_exposer.def( "__len__", &__len_count< ::SireFF::EnergyTable > );
    }

}
