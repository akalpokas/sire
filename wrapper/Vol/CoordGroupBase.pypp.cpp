// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CoordGroupBase.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "coordgroup.h"

#include <QDebug>

#include "coordgroup.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_CoordGroupBase_class(){

    { //::SireVol::CoordGroupBase
        typedef bp::class_< SireVol::CoordGroupBase, boost::noncopyable > CoordGroupBase_exposer_t;
        CoordGroupBase_exposer_t CoordGroupBase_exposer = CoordGroupBase_exposer_t( "CoordGroupBase", "This is the base class of all CoordGroup-like classes\n(e.g. CoordGroup and CoordGroupEditor). CoordGroup classes\nhold a group of coordinates, together with an axis-aligned\nbox that completely encloses all of those points. The\nclass is implicitly shared, and, since it is used in the\nmost performance-sensitive parts of the code, has\na slightly more complex implicit-sharing design.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope CoordGroupBase_scope( CoordGroupBase_exposer );
        { //::SireVol::CoordGroupBase::aaBox
        
            typedef ::SireVol::AABox const & ( ::SireVol::CoordGroupBase::*aaBox_function_type)(  ) const;
            aaBox_function_type aaBox_function_value( &::SireVol::CoordGroupBase::aaBox );
            
            CoordGroupBase_exposer.def( 
                "aaBox"
                , aaBox_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the enclosing AABox" );
        
        }
        { //::SireVol::CoordGroupBase::assertSameSize
        
            typedef void ( ::SireVol::CoordGroupBase::*assertSameSize_function_type)( ::QVector< SireMaths::Vector > const & ) const;
            assertSameSize_function_type assertSameSize_function_value( &::SireVol::CoordGroupBase::assertSameSize );
            
            CoordGroupBase_exposer.def( 
                "assertSameSize"
                , assertSameSize_function_value
                , ( bp::arg("coordinates") )
                , bp::release_gil_policy()
                , "Assert that the array of coordinates contains the same number of\ncoordinates as are in this group\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupBase::assertSameSize
        
            typedef void ( ::SireVol::CoordGroupBase::*assertSameSize_function_type)( ::SireVol::CoordGroupBase const & ) const;
            assertSameSize_function_type assertSameSize_function_value( &::SireVol::CoordGroupBase::assertSameSize );
            
            CoordGroupBase_exposer.def( 
                "assertSameSize"
                , assertSameSize_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Assert that the CoordGroup other contains the same number of\ncoordinates as are in this group\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupBase::assertValidIndex
        
            typedef void ( ::SireVol::CoordGroupBase::*assertValidIndex_function_type)( ::quint32 ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireVol::CoordGroupBase::assertValidIndex );
            
            CoordGroupBase_exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Assert that i is a valid index in this group\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupBase::at
        
            typedef ::SireMaths::Vector const & ( ::SireVol::CoordGroupBase::*at_function_type)( ::quint32 ) const;
            at_function_type at_function_value( &::SireVol::CoordGroupBase::at );
            
            CoordGroupBase_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the ith coordinate of this group\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupBase::count
        
            typedef int ( ::SireVol::CoordGroupBase::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireVol::CoordGroupBase::count );
            
            CoordGroupBase_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of coordinates in this group" );
        
        }
        { //::SireVol::CoordGroupBase::isEmpty
        
            typedef bool ( ::SireVol::CoordGroupBase::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireVol::CoordGroupBase::isEmpty );
            
            CoordGroupBase_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether this group is empty (has no coordinates)" );
        
        }
        { //::SireVol::CoordGroupBase::maybeDifferent
        
            typedef bool ( ::SireVol::CoordGroupBase::*maybeDifferent_function_type)( ::SireVol::CoordGroupBase const & ) const;
            maybeDifferent_function_type maybeDifferent_function_value( &::SireVol::CoordGroupBase::maybeDifferent );
            
            CoordGroupBase_exposer.def( 
                "maybeDifferent"
                , maybeDifferent_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return whether other may be different to this group - this uses\na simple comparison of the memory addresses of the storage of\nthese two groups to see if they are definitely the same, or maybe\ndifferent" );
        
        }
        CoordGroupBase_exposer.def( bp::self != bp::self );
        CoordGroupBase_exposer.def( bp::self == bp::self );
        { //::SireVol::CoordGroupBase::operator[]
        
            typedef ::SireMaths::Vector const & ( ::SireVol::CoordGroupBase::*__getitem___function_type)( ::quint32 ) const;
            __getitem___function_type __getitem___function_value( &::SireVol::CoordGroupBase::operator[] );
            
            CoordGroupBase_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireVol::CoordGroupBase::size
        
            typedef int ( ::SireVol::CoordGroupBase::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireVol::CoordGroupBase::size );
            
            CoordGroupBase_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of coordinates in this group" );
        
        }
        { //::SireVol::CoordGroupBase::toString
        
            typedef ::QString ( ::SireVol::CoordGroupBase::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireVol::CoordGroupBase::toString );
            
            CoordGroupBase_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this CoordGroup" );
        
        }
        { //::SireVol::CoordGroupBase::toVector
        
            typedef ::QVector< SireMaths::Vector > ( ::SireVol::CoordGroupBase::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireVol::CoordGroupBase::toVector );
            
            CoordGroupBase_exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "Return an array containing just the coordinates in this group" );
        
        }
        { //::SireVol::CoordGroupBase::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::CoordGroupBase::typeName );
            
            CoordGroupBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireVol::CoordGroupBase::what
        
            typedef char const * ( ::SireVol::CoordGroupBase::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireVol::CoordGroupBase::what );
            
            CoordGroupBase_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CoordGroupBase_exposer.staticmethod( "typeName" );
        CoordGroupBase_exposer.def( "__str__", &__str__< ::SireVol::CoordGroupBase > );
        CoordGroupBase_exposer.def( "__repr__", &__str__< ::SireVol::CoordGroupBase > );
        CoordGroupBase_exposer.def( "__len__", &__len_size< ::SireVol::CoordGroupBase > );
    }

}
