--09-28-05 At's tf information

create schema at_tfdb;

set search_path to at_tfdb;

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
	comment	varchar,
	source	varchar
	);

create table family(
	id	serial primary key,
	family_id	varchar unique,
	description	varchar
	);

create table family_ref(
	id	serial primary key,
	family_id	varchar references family(family_id) on delete cascade on update cascade,
	ref_title       varchar,
	ref_authors     varchar,
	ref_external_link       varchar
        );
	
create table factor(
	id	serial primary key,
	tf_acc	varchar unique,
	tf_id	varchar,
	tf_name	varchar,
	tf_syn	varchar,
	organism	varchar,
	peptide	varchar,
	sequence	varchar,
	sequence_source	varchar,
	families	varchar[],	--list of family_id from table family
	comment	varchar
	);

create table factor_seq(
	tf_acc	varchar,
	sequence	varchar, --DNA
	peptide	varchar
	);

create table matrix(
	id	serial primary key,
	mt_acc	varchar unique,
	mt_id	varchar unique,
	family_id	varchar,
	pwm	float[][],
	consensus	varchar,
	comment	varchar
	);


create table binding_site(
	id	serial primary key,
	bs_id	varchar	unique,
	mt_id	varchar references matrix(mt_id) on delete cascade on update cascade,
	chromosome	varchar,
	strand	varchar(1),
	bs_genome_start	integer,
	bs_genome_end	integer,
	bs_disp_start	integer,
	bs_disp_end	integer,
	prom_id	integer references prom_seq(id) on delete cascade on update cascade,
	color	varchar,
	matrix_similarity_score float,
	sequence	varchar,
	motif 	varchar
	);

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

--09-28-05 different from transfacdb.sql's same function
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
		
		IF NEW.bs_genome_start IS NOT NULL AND NEW.bs_genome_end IS NOT NULL AND NEW.sequence IS NOT NULL THEN
			IF NEW.bs_disp_start IS NULL OR NEW.bs_disp_end IS NULL THEN
				select into prom_seq_position * from find_prom_seq_position(NEW.prom_id)
				as (prom_genome_start integer, chromosome varchar);
				IF prom_seq_position.prom_genome_start IS NOT NULL THEN
					NEW.bs_disp_start := NEW.bs_genome_start - prom_seq_position.prom_genome_start +1;
					NEW.bs_disp_end := NEW.bs_disp_start + length(NEW.sequence) - 1;
				END IF;
			END IF;
		END IF;
		RETURN NEW;
	END;
	' LANGUAGE plpgsql;

CREATE TRIGGER before_binding_site BEFORE INSERT OR UPDATE ON binding_site
	FOR EACH ROW EXECUTE PROCEDURE before_binding_site();

