use strict;
use warnings;


my @files=glob("*.tab");

foreach my $file (@files){
    print STDERR $file, "\n";
    my $newName=$file;
    $newName=~s/Lister2009/Listeretal_Nature_2009/;
    $newName=~s/Lister2011/Listeretal_Nature_2011/;
    $newName=~s/Hannon2011/Hodgesetal_MolCell_2011/;
    print STDERR $newName, "\n";
    print "\n";
    system("mv $file $newName");
}
