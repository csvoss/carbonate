# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Property'
        db.create_table(u'api_property', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50)),
        ))
        db.send_create_signal(u'api', ['Property'])

        # Adding model 'Reagent'
        db.create_table(u'api_reagent', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('is_solvent', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('diagram_name', self.gf('django.db.models.fields.CharField')(max_length=50, blank=True)),
            ('smiles', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
        ))
        db.send_create_signal(u'api', ['Reagent'])

        # Adding M2M table for field properties on 'Reagent'
        m2m_table_name = db.shorten_name(u'api_reagent_properties')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reagent', models.ForeignKey(orm[u'api.reagent'], null=False)),
            ('property', models.ForeignKey(orm[u'api.property'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reagent_id', 'property_id'])

        # Adding model 'ReagentSet'
        db.create_table(u'api_reagentset', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=50, blank=True)),
            ('solvent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['api.Reagent'], null=True, blank=True)),
        ))
        db.send_create_signal(u'api', ['ReagentSet'])

        # Adding M2M table for field reagents on 'ReagentSet'
        m2m_table_name = db.shorten_name(u'api_reagentset_reagents')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reagentset', models.ForeignKey(orm[u'api.reagentset'], null=False)),
            ('reagent', models.ForeignKey(orm[u'api.reagent'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reagentset_id', 'reagent_id'])

        # Adding M2M table for field solvent_properties on 'ReagentSet'
        m2m_table_name = db.shorten_name(u'api_reagentset_solvent_properties')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('reagentset', models.ForeignKey(orm[u'api.reagentset'], null=False)),
            ('property', models.ForeignKey(orm[u'api.property'], null=False))
        ))
        db.create_unique(m2m_table_name, ['reagentset_id', 'property_id'])

        # Adding model 'Reaction'
        db.create_table(u'api_reaction', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('process_function', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('reagent_set', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactions', to=orm['api.ReagentSet'])),
        ))
        db.send_create_signal(u'api', ['Reaction'])


    def backwards(self, orm):
        # Deleting model 'Property'
        db.delete_table(u'api_property')

        # Deleting model 'Reagent'
        db.delete_table(u'api_reagent')

        # Removing M2M table for field properties on 'Reagent'
        db.delete_table(db.shorten_name(u'api_reagent_properties'))

        # Deleting model 'ReagentSet'
        db.delete_table(u'api_reagentset')

        # Removing M2M table for field reagents on 'ReagentSet'
        db.delete_table(db.shorten_name(u'api_reagentset_reagents'))

        # Removing M2M table for field solvent_properties on 'ReagentSet'
        db.delete_table(db.shorten_name(u'api_reagentset_solvent_properties'))

        # Deleting model 'Reaction'
        db.delete_table(u'api_reaction')


    models = {
        u'api.property': {
            'Meta': {'object_name': 'Property'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'api.reaction': {
            'Meta': {'object_name': 'Reaction'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'process_function': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'reagent_set': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactions'", 'to': u"orm['api.ReagentSet']"})
        },
        u'api.reagent': {
            'Meta': {'object_name': 'Reagent'},
            'diagram_name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_solvent': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'properties': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['api.Property']", 'null': 'True', 'blank': 'True'}),
            'smiles': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'})
        },
        u'api.reagentset': {
            'Meta': {'object_name': 'ReagentSet'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'blank': 'True'}),
            'reagents': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'reagentSets'", 'symmetrical': 'False', 'to': u"orm['api.Reagent']"}),
            'solvent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['api.Reagent']", 'null': 'True', 'blank': 'True'}),
            'solvent_properties': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['api.Property']", 'null': 'True', 'blank': 'True'})
        }
    }

    complete_apps = ['api']