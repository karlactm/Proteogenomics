package FASTAParser;
use Protein;

sub new {
    my $class = shift;
    my $self = {
        _filePath => shift,
        _proteins => []
    };
    bless $self, $class;
    return $self;
}

sub parse {
    my $self = shift;
    my $path = $self->{_filePath};
    open (my $file, "<", $path);
    my $id = "";
    my $description = "";
    my $sample = "";
    my $header = "";
    my $aminoacids = "";
    my $i = 0;
    while (my $line = <$file>) {
	chomp ($line);
	if ($line =~ m/^>/) {
        $i++;
	    if ($header ne ""){
	        my $fasta = Protein->new($id, $description, $sample, $header, $aminoacids);
		    push @{$self->{_proteins}}, $fasta;
		    $aminoacids = "";	
            $description = "";
            $sample = "";
	     }
	     my @array = split(">",$line);
	     $header = $array[1]; 
             @array = split(/\s/,$header);
             $id = $array[0];
             for (my $i = 1; $i < scalar (@array); $i++) {
                 $description .= $array[$i]." ";
             }
             @array = split("", $header);
             my $position = 0;
             for (my $i = 0; $i < scalar (@array); $i++) {
                 if ($array [$i] eq "[") {
                     $position = $i;
                 }
             }
             for (my $i = ($position + 1); $i < (scalar (@array) - 1); $i++) {
                 $sample .= $array [$i];
             }
        }
        else{
            $aminoacids.= $line;
        }
    }
    my $fasta = Protein->new($id, $description, $sample, $header, $aminoacids);
    push @{$self->{_proteins}}, $fasta;
    return 1;	
}

sub getProteins {
    my $self = shift;
    return \@{$self->{_proteins}};
}

sub getProteinsSize {
    my $self = shift;
    return scalar @{$self->{_proteins}};
}

sub getFilePath {
    my $self = shift;
    return $self->{_filePath};
}

sub getFileName {
    my $self = shift;
    my $path = $self->{_filePath};
    my @fileName = ();
    foreach (split(/\//,$path)) {
        if ($_ =~ m/^[GCA]/gci) { 
            @fileName = split(/\_/,$_);
        }
    } 
    return $fileName[1];  
}

1;
