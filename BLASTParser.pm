package BLASTParser;
use BLASTProtein;

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
    open(my $file, "<", $path) or die "Não foi possível abrir o arquivo '$path' $!";
    my $query_id = undef;
    my $subject_id = undef;
    my $identity = undef;
    my $positives = undef;
    my $bit_score = undef;
    my $subject_cov = undef;
    my $hsp_cov = undef;
    my $count = 0;
    my @query = ();
    while (my $line = <$file>) {
	chomp ($line);
	if($line =~ m/^#/) {
	    $count++;
	    if ($count == 2) {
	        @query = split (/\s/, $line);
	    }
	    if ($line =~ /^# 0 hits found$/) {
	           my $blast = BLASTProtein->new($query[2], undef, 0, 0, 0, 0, 0);
	           push @{$self->{_proteins}}, $blast;
	           $count = 0;
	           @query = ();
	    }	
	} else {
             $count = 0;	
	     my @fields = split (/\t/,$line);
             $query_id = $fields[0];
             $subject_id = $fields[1];
             $identity = $fields[2];
             $positives = $fields[3];
             $bit_score = $fields[12];
             $subject_cov = $fields[15];
             $hsp_cov = $fields[16];
             my $blast = BLASTProtein->new($query_id, $subject_id, $identity, $positives, $bit_score, $subject_cov, $hsp_cov);
             push @{$self->{_proteins}}, $blast;
       } 
     }
    return 1;	
}

sub getProteins {
    my $self = shift;
    return \@{$self->{_proteins}};
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
        if ($_ =~ m/[0-9]/) { 
            $fileName = $_;
        }
    } 
    return $fileName;  
}

1;
