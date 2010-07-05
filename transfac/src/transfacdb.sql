create schema transfac;

set search_path to transfac;

create table readme(
	id      serial primary key,
	name    varchar,
	description     varchar,
	date_created	timestamp default current_timestamp
);
				

create table prom_type(
	id	integer primary key,
	name	varchar
	);

insert into prom_type values(1, 'upstream');
insert into prom_type values(2, '5UTR');
insert into prom_type values(3, '3UTR');
insert into prom_type values(4, 'downstream');

create table prom_seq(
	id	serial primary key,
	prom_acc	varchar not null,
	chromosome	varchar,
	strand	varchar(1),
	prom_genome_start	integer,
	prom_genome_end	integer,
	organism	varchar,
	prom_type_id	integer references prom_type(id) on delete cascade on update cascade,
	sequence	varchar not null,
	masked_seq	varchar,
	comment	varchar
	);

--01-03-06 add index
create index prom_seq_prom_acc_idx on prom_seq(prom_acc);

--function arguments($1,$2...) are constant, can't be modified.
create function insert_prom_seq(varchar, varchar, varchar, integer, integer,
	varchar, integer, varchar, varchar) returns integer as '
	declare
		_prom_acc varchar;
		_chromosome varchar;
		_strand varchar(1);
		_prom_genome_start integer;
		_prom_genome_end integer;
		_organism varchar;
		_prom_type_id integer;
		_sequence varchar;
		_comment varchar;
	begin
		_prom_acc := $1;
		_chromosome := $2;
		_strand := $3;
		_prom_genome_start := $4;
		_prom_genome_end := $5;
		_organism := $6;
		_prom_type_id := $7;
		_sequence := $8;
		_comment := $9;
		if _prom_genome_start > _prom_genome_end then
			return 1;
		end if;
		if _prom_acc = '''' then
			_prom_acc := NULL;
		end if;
		if _chromosome ='''' then
			_chromosome := NULL;
		end if;
		if _strand = '''' then
			_strand := NULL;
		end if;
		if _organism = '''' then
			_organism := NULL;
		end if;
		if _sequence = '''' then
			_sequence := NULL;
		end if;
		if _comment = '''' then
			_comment := NULL;
		end if;
		insert into prom_seq(prom_acc, chromosome, strand, prom_genome_start, prom_genome_end, organism, prom_type_id, sequence, comment)
			values(_prom_acc, _chromosome, _strand, _prom_genome_start, _prom_genome_end, _organism, _prom_type_id, _sequence, _comment);
		if FOUND then
			return 0;
		else
			return 2;
		end if;
	end;
	'
	language plpgsql;

create table reference(
	id	serial primary key,
	ref_acc	varchar unique,
	ref_external_link	varchar,
	ref_authors	varchar,
	ref_title	varchar,
	ref_journal	varchar
	);

CREATE FUNCTION before_reference() RETURNS trigger AS '
    BEGIN
        IF NEW.ref_acc='''' THEN
            NEW.ref_acc := null;
        END IF;
        IF NEW.ref_external_link='''' THEN
            NEW.ref_external_link := null;
        END IF;
	IF NEW.ref_authors='''' THEN
		NEW.ref_authors := null;
	END IF;
	IF NEW.ref_title='''' THEN
		NEW.ref_title := null;
	END IF;
	IF NEW.ref_journal = ''{}'' THEN
		NEW.ref_journal := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_reference BEFORE INSERT OR UPDATE ON reference
	FOR EACH ROW EXECUTE PROCEDURE before_reference();

-- 2009-9-9 store the new families_data.tbl downloaded from the AGRIS site --
create table families_data20090909 (
FamilyID 	varchar,
LocusName 	varchar,
GeneName  varchar,
Description 	varchar,
RNIntegration 	varchar,
SubFamily	 varchar,
BindingSite	 varchar,
SpecialInfo	 varchar,
unknownField  varchar
);


	
create table factor(
	id	serial primary key,
	tf_acc	varchar unique,
	tf_id	varchar,
	tf_name	varchar,
	tf_syn	varchar,
	organism	varchar,
	gene_acc	varchar,
	homolog	varchar,
	class_acc	varchar,
	class_tree_id	varchar,
	size	varchar,
	sequence	varchar,
	sequence_source	varchar,
	feature	varchar[],
	struc_feature	varchar[],
	cell_positive	varchar,
	cell_negative	varchar,
	expr_pattern	varchar[],
	func_feature	varchar[],
	interacting_partners	varchar[],
	matrices	varchar[],
	sites	varchar[],
	binding_region	varchar[],
	external_database_links	varchar[],
	reference	varchar[]
	);


CREATE OR REPLACE FUNCTION before_factor() RETURNS trigger AS '
    BEGIN
        IF NEW.tf_acc='''' THEN
            NEW.tf_acc := null;
        END IF;
        IF NEW.tf_id='''' THEN
            NEW.tf_id := null;
        END IF;
	IF NEW.tf_name='''' THEN
		NEW.tf_name := null;
	END IF;
	IF NEW.tf_syn='''' THEN
		NEW.tf_syn := null;
	END IF;
	IF NEW.organism = '''' THEN
		NEW.organism := null;
	END IF;
	IF NEW.gene_acc = '''' THEN
		NEW.gene_acc := null;
	END IF;
	IF NEW.homolog = '''' THEN
		NEW.homolog := null;
	END IF;
	IF NEW.class_acc = '''' THEN
		NEW.class_acc := null;
	END IF;
	IF NEW.class_tree_id = '''' THEN
		NEW.class_tree_id := null;
	END IF;
	IF NEW.size = '''' THEN
		NEW.size := null;
	END IF;	
	IF NEW.sequence = '''' THEN
		NEW.sequence := null;
	END IF;
	IF NEW.sequence_source = '''' THEN
		NEW.sequence_source := null;
	END IF;
	IF NEW.feature = ''{}'' THEN
		NEW.feature := null;
	END IF;
	IF NEW.struc_feature = ''{}'' THEN
		NEW.struc_feature := null;
	END IF;
	IF NEW.cell_positive = '''' THEN
		NEW.cell_positive := null;
	END IF;
	IF NEW.cell_negative = '''' THEN
		NEW.cell_negative := null;
	END IF;
	IF NEW.expr_pattern = ''{}'' THEN
		NEW.expr_pattern := null;
	END IF;
	IF NEW.func_feature = ''{}'' THEN
		NEW.func_feature := null;
	END IF;
	IF NEW.interacting_partners = ''{}'' THEN
		NEW.interacting_partners := null;
	END IF;
	IF NEW.matrices = ''{}'' THEN
		NEW.matrices := null;
	END IF;
	IF NEW.sites = ''{}'' THEN
		NEW.sites := null;
	END IF;
	IF NEW.binding_region = ''{}'' THEN
		NEW.binding_region := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
        IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_factor BEFORE INSERT OR UPDATE ON factor
	FOR EACH ROW EXECUTE PROCEDURE before_factor();


create table matrix(
	id	serial primary key,
	mt_acc	varchar unique,
	mt_id	varchar unique,
	tf_name	varchar,
	tf_description	varchar,
	factors	varchar[],
	pwm	float[][],
	consensus	varchar,
	basis	varchar,
	sites	varchar[],
	site_in_matrix_accs	varchar[],
	comment	varchar,
	reference	varchar[]
	);

create table site_in_matrix(
	id	serial primary key,
	acc	varchar unique,
	mt_acc	varchar,
	site_acc	varchar,
	sequence	varchar,
	mt_start	integer,
	len	integer,
	gap_list	integer[],
	orientation	varchar(1)
	);

create view matrix2factor as select m.mt_acc, m.mt_id, m.tf_name as factor_name, m.tf_description, m.pwm, 
	m.consensus, m.basis, m.sites as mt_sites,m.comment as mt_comment, m.reference as mt_reference, 
	f.tf_acc, f.tf_id, f.tf_name, f.tf_syn, f.organism, f.gene_acc, f.homolog, f.class_acc, f.class_tree_id,
	f.size, f.sequence, f.sequence_source, f.feature, f.struc_feature, f.cell_positive, f.cell_negative,
	f.expr_pattern, f.func_feature, f.interacting_partners, f.sites as tf_sites, f.binding_region, f.external_database_links,
	f.reference as tf_reference from matrix m, factor f
	where f.tf_acc=any(m.factors);

create view factor2matrix as select f.tf_acc, f.tf_id, m.mt_acc, m.mt_id from matrix m, factor f
	where m.mt_acc=any(f.matrices);


CREATE OR REPLACE FUNCTION before_matrix() RETURNS trigger AS '
    BEGIN
        IF NEW.mt_acc='''' THEN
            NEW.mt_acc := null;
        END IF;
        IF NEW.mt_id='''' THEN
            NEW.mt_id := null;
        END IF;
	IF NEW.tf_name='''' THEN
		NEW.tf_name := null;
	END IF;
	IF NEW.tf_description='''' THEN
		NEW.tf_description := null;
	END IF;
	IF NEW.factors = ''{}'' THEN
		NEW.factors := null;
	END IF;
	IF NEW.pwm = ''{}'' THEN
		NEW.pwm := null;
	END IF;
	IF NEW.consensus = '''' THEN
		NEW.consensus := null;
	END IF;
	IF NEW.basis = '''' THEN
		NEW.basis := null;
	END IF;
	IF NEW.sites = ''{}'' THEN
		NEW.sites := null;
	END IF;
	IF NEW.site_in_matrix_accs = ''{}'' THEN
		NEW.site_in_matrix_accs := null;
	END IF;
	IF NEW.comment = '''' THEN
		NEW.comment := null;
	END IF;
	IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_matrix BEFORE INSERT OR UPDATE ON matrix
	FOR EACH ROW EXECUTE PROCEDURE before_matrix();


CREATE OR REPLACE FUNCTION before_site_in_matrix() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
	IF NEW.mt_acc='''' THEN
            NEW.mt_acc := null;
        END IF;
        IF NEW.site_acc='''' THEN
            NEW.site_acc := null;
        END IF;
	IF NEW.sequence = '''' THEN
		NEW.sequence := null;
	END IF;
	IF NEW.mt_start=2000000 THEN
		NEW.mt_start := null;
	END IF;
	IF NEW.len=2000000 THEN
		NEW.len := null;
	END IF;
	IF NEW.gap_list = ''{}'' THEN
		NEW.gap_list := null;
	END IF;
	IF NEW.orientation = '''' THEN
		NEW.orientation := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_site_in_matrix BEFORE INSERT OR UPDATE ON site_in_matrix
	FOR EACH ROW EXECUTE PROCEDURE before_site_in_matrix();

create table binding_site(
	id	serial primary key,
	mt_id	varchar references matrix(mt_id) on delete cascade on update cascade,
	chromosome	varchar,
	strand	varchar(1),
	bs_genome_start	integer,
	bs_genome_end	integer,
	bs_disp_start	integer,
	bs_disp_end	integer,
	prom_id	integer references prom_seq(id) on delete cascade on update cascade,
	core_similarity_score	float,
	matrix_similarity_score	float,
	sequence	varchar
	);

--01-03-06 add index
create INDEX binding_site_prom_id_idx on binding_site(prom_id);

--used by before_binding_site()
CREATE OR REPLACE FUNCTION find_prom_seq_position(integer) RETURNS record AS '
	DECLARE
		_prom_id ALIAS FOR $1;
		result RECORD;
	BEGIN
		select prom_genome_start, chromosome into result from prom_seq where  id=_prom_id;
		RETURN result;
	END;
	'
	LANGUAGE plpgsql
	STABLE
	RETURNS NULL ON NULL INPUT;
	
CREATE OR REPLACE FUNCTION before_binding_site() RETURNS trigger AS '
	DECLARE
		prom_seq_position RECORD;
	BEGIN
		IF NEW.mt_id='''' THEN
			NEW.mt_id := null;
		END IF;
		IF NEW.chromosome IS NULL THEN
			select into prom_seq_position * from find_prom_seq_position(NEW.prom_id)
				as (prom_genome_start integer, chromosome varchar);
			IF prom_seq_position.chromosome IS NOT NULL THEN
				NEW.chromosome := prom_seq_position.chromosome;
			ELSE
				RETURN null;	--abort this action
			END IF;
		END IF;
		IF NEW.strand='''' THEN
			NEW.strand := null;
		END IF;
		IF NEW.bs_disp_start < 0 THEN
			NEW.bs_disp_start := null;
		END IF;
		IF NEW.bs_disp_end < 0 THEN
			NEW.bs_disp_end := null;
		END IF;
		IF NEW.sequence = '''' THEN
			NEW.sequence := null;
		END IF;
		
		IF NEW.bs_genome_start IS NULL or NEW.bs_genome_end IS NULL THEN
			IF NEW.bs_disp_start IS NOT NULL and NEW.sequence IS NOT NULL THEN
				select into prom_seq_position * from find_prom_seq_position(NEW.prom_id)
				as (prom_genome_start integer, chromosome varchar);
				IF prom_seq_position.prom_genome_start IS NOT NULL THEN
					NEW.bs_genome_start := prom_seq_position.prom_genome_start + NEW.bs_disp_start -1;
					NEW.bs_genome_end := NEW.bs_genome_start + length(NEW.sequence) - 1;
				END IF;
			END IF;
		END IF;
		RETURN NEW;
	END;
	' LANGUAGE plpgsql;

CREATE TRIGGER before_binding_site BEFORE INSERT OR UPDATE ON binding_site
	FOR EACH ROW EXECUTE PROCEDURE before_binding_site();

CREATE TABLE site(
	id	serial primary key,
	acc	varchar unique,
	acc_id	varchar,
	type	varchar,
	description	varchar,
	gene_acc	varchar,
	organism	varchar,
	gene_region	varchar,
	sequence	varchar,
	element	varchar,
	position_ref_point	varchar,
	position_start	integer,
	position_end	integer,
	cell_source	varchar[],
	method	varchar[],
	comment	varchar,
	external_database_links	varchar[],
	reference	varchar[]
	);

create view matrix2site as select m.mt_acc, m.mt_id, m.tf_name as factor_name, m.tf_description, m.pwm, 
	m.consensus, m.basis, m.sites as mt_sites,m.comment as mt_comment, m.reference as mt_reference,
	s.acc as site_acc, s.acc_id as site_id, s.type, s.description, s.gene_acc, s.organism, s.gene_region, 
	s.sequence, s.element, s.position_ref_point, s.position_start, s.position_end, s.cell_source, s.method,
	s.comment as site_comment, s.external_database_links, s.reference as site_reference, sm.acc,
	sm.sequence as site_in_matrix_seq, sm.mt_start, sm.len as mt_len, sm.gap_list, sm.orientation
	from matrix m, site s, site_in_matrix sm
	where s.acc=any(m.sites) and sm.mt_acc=m.mt_acc and sm.site_acc=s.acc;

CREATE OR REPLACE FUNCTION before_site() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
        IF NEW.acc_id='''' THEN
            NEW.acc_id := null;
        END IF;
	IF NEW.type ='''' THEN
		NEW.type := null;
	END IF;
	IF NEW.description='''' THEN
		NEW.description := null;
	END IF;
	IF NEW.gene_acc='''' THEN
		NEW.gene_acc := null;
	END IF;
	IF NEW.organism='''' THEN
		NEW.organism := null;
	END IF;
	IF NEW.gene_region='''' THEN
		NEW.gene_region := null;
	END IF;
	IF NEW.sequence='''' THEN
		NEW.sequence := null;
	END IF;
	IF NEW.element='''' THEN
		NEW.element:= null;
	END IF;
	IF NEW.position_ref_point='''' THEN
		NEW.position_ref_point := null;
	END IF;
	IF NEW.position_start=2000000 THEN
		NEW.position_start := null;
	END IF;
	IF NEW.position_end=2000000 THEN
		NEW.position_end := null;
	END IF;
	
	IF NEW.cell_source = ''{}'' THEN
		NEW.cell_source := null;
	END IF;
	IF NEW.method = ''{}'' THEN
		NEW.method := null;
	END IF;
	IF NEW.comment = '''' THEN
		NEW.comment := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
	IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_site BEFORE INSERT OR UPDATE ON site
	FOR EACH ROW EXECUTE PROCEDURE before_site();

CREATE TABLE gene(
	id	serial primary key,
	acc	varchar unique,
	acc_id	varchar,
	short_description	varchar,
	description	varchar,
	synonyms	varchar,
	organism	varchar,
	chromosome	varchar,
	prom_classification	varchar,
	regulation	varchar,
	binding_region	varchar[],
	factors	varchar[],
	external_database_links	varchar[],
	reference	varchar[]
	);

CREATE OR REPLACE FUNCTION before_gene() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
        IF NEW.acc_id='''' THEN
            NEW.acc_id := null;
        END IF;
	IF NEW.short_description ='''' THEN
		NEW.short_description := null;
	END IF;
	IF NEW.description='''' THEN
		NEW.description := null;
	END IF;
	IF NEW.synonyms='''' THEN
		NEW.synonyms := null;
	END IF;
	IF NEW.organism='''' THEN
		NEW.organism := null;
	END IF;
	IF NEW.chromosome='''' THEN
		NEW.chromosome := null;
	END IF;
	IF NEW.prom_classification='''' THEN
		NEW.prom_classification := null;
	END IF;
	IF NEW.regulation ='''' THEN
		NEW.regulation := null;
	END IF;
	IF NEW.binding_region=''{}'' THEN
		NEW.binding_region := null;
	END IF;
	IF NEW.factors = ''{}'' THEN
		NEW.factors := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
	IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_gene BEFORE INSERT OR UPDATE ON gene
	FOR EACH ROW EXECUTE PROCEDURE before_gene();

CREATE TABLE cell(
	id	serial primary key,
	acc	varchar unique,
	cell_factor_source	varchar,
	organism	varchar,
	cell_description varchar,
	sites	varchar[],
	binding_region	varchar[],
	external_database_links	varchar[]
	);

CREATE OR REPLACE FUNCTION before_cell() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
        IF NEW.cell_factor_source='''' THEN
            NEW.cell_factor_source := null;
        END IF;
	IF NEW.organism='''' THEN
		NEW.organism := null;
	END IF;
	IF NEW.cell_description='''' THEN
		NEW.cell_description := null;
	END IF;
	IF NEW.sites=''{}'' THEN
		NEW.sites := null;
	END IF;
	IF NEW.binding_region =''{}'' THEN
		NEW.binding_region := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_cell BEFORE INSERT OR UPDATE ON cell
	FOR EACH ROW EXECUTE PROCEDURE before_cell();

CREATE TABLE fragment(
	id	serial primary key,
	acc	varchar,
	acc_id	varchar,
	description	varchar,
	gene_acc	varchar,
	organism	varchar,
	sequence	varchar,
	sequence_source	varchar,
	factors	varchar[],
	method	varchar[],
	external_database_links	varchar[],
	reference	varchar[]
	);

CREATE OR REPLACE FUNCTION before_fragment() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
	 IF NEW.acc_id='''' THEN
            NEW.acc_id := null;
        END IF;
	IF NEW.description='''' THEN
		NEW.description := null;
	END IF;
	IF NEW.gene_acc='''' THEN
		NEW.gene_acc := null;
	END IF;
	IF NEW.organism='''' THEN
		NEW.organism := null;
	END IF;
	IF NEW.sequence='''' THEN
		NEW.sequence := null;
	END IF;
	IF NEW.sequence_source='''' THEN
		NEW.sequence_source := null;
	END IF;
	IF NEW.factors=''{}'' THEN
		NEW.factors := null;
	END IF;
	IF NEW.method =''{}'' THEN
		NEW.method := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
	IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_fragment BEFORE INSERT OR UPDATE ON fragment
	FOR EACH ROW EXECUTE PROCEDURE before_fragment();

CREATE TABLE class(
	id	serial primary key,
	acc	varchar,
	acc_id	varchar,
	class_struc	varchar,
	description	varchar,
	factors	varchar[],
	external_database_links	varchar[],
	reference	varchar[]
	);

CREATE OR REPLACE FUNCTION before_class() RETURNS trigger AS '
    BEGIN
        IF NEW.acc='''' THEN
            NEW.acc := null;
        END IF;
	 IF NEW.acc_id='''' THEN
            NEW.acc_id := null;
        END IF;
	IF NEW.class_struc='''' THEN
		NEW.class_struc := null;
	END IF;
	IF NEW.description='''' THEN
		NEW.description := null;
	END IF;
	IF NEW.factors=''{}'' THEN
		NEW.factors := null;
	END IF;
	IF NEW.external_database_links = ''{}'' THEN
		NEW.external_database_links := null;
	END IF;
	IF NEW.reference = ''{}'' THEN
        	NEW.reference := null;
	END IF;
        RETURN NEW;
    END;
' LANGUAGE plpgsql;

CREATE TRIGGER before_class BEFORE INSERT OR UPDATE ON class
	FOR EACH ROW EXECUTE PROCEDURE before_class();

--01-16-06
create table matrix2no_of_random_hits(
	id	serial primary key,
	mt_id	varchar,
	gc_perc	float,
	no_of_hits	integer,
	seq_id	integer
	);

--01-31-06 LinkFactorID2EntrezGeneID.py fills in these two tables
create table tf_acc2entrezgene_id_raw(
	id	serial primary key,
	tf_acc	varchar,
	gene_id	integer,
	bridge_acc	varchar
	);

create table tf_acc2entrezgene_id(
	id	serial primary key,
	tf_acc	varchar,
	gene_id	integer
	);

CREATE OR REPLACE VIEW mt2gene_id AS SELECT distinct m.mt_acc, m.mt_id, t.gene_id, g.tax_id from matrix m, tf_acc2entrezgene_id t, gene.gene g where t.tf_acc = any(m.factors) and g.gene_id=t.gene_id;
