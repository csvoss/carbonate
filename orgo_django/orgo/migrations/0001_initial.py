# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'MoleculeBoxModel'
        db.create_table(u'orgo_moleculeboxmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('moleculeBox', self.gf('orgo.models.PickledObjectField')(null=True)),
            ('svg', self.gf('django.db.models.fields.TextField')(null=True)),
            ('equalsTarget', self.gf('django.db.models.fields.BooleanField')()),
            ('isStartingMaterial', self.gf('django.db.models.fields.BooleanField')()),
        ))
        db.send_create_signal(u'orgo', ['MoleculeBoxModel'])

        # Adding model 'ArrowModel'
        db.create_table(u'orgo_arrowmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('pointFrom', self.gf('django.db.models.fields.related.ForeignKey')(related_name='arrowPointsFrom', null=True, to=orm['orgo.MoleculeBoxModel'])),
            ('pointTo', self.gf('django.db.models.fields.related.ForeignKey')(related_name='arrowPointsTo', null=True, to=orm['orgo.MoleculeBoxModel'])),
            ('reagentsHtml', self.gf('django.db.models.fields.TextField')(null=True)),
        ))
        db.send_create_signal(u'orgo', ['ArrowModel'])

        # Adding model 'SolutionModel'
        db.create_table(u'orgo_solutionmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
        ))
        db.send_create_signal(u'orgo', ['SolutionModel'])

        # Adding M2M table for field molecules on 'SolutionModel'
        m2m_table_name = db.shorten_name(u'orgo_solutionmodel_molecules')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('solutionmodel', models.ForeignKey(orm[u'orgo.solutionmodel'], null=False)),
            ('moleculeboxmodel', models.ForeignKey(orm[u'orgo.moleculeboxmodel'], null=False))
        ))
        db.create_unique(m2m_table_name, ['solutionmodel_id', 'moleculeboxmodel_id'])

        # Adding M2M table for field arrows on 'SolutionModel'
        m2m_table_name = db.shorten_name(u'orgo_solutionmodel_arrows')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('solutionmodel', models.ForeignKey(orm[u'orgo.solutionmodel'], null=False)),
            ('arrowmodel', models.ForeignKey(orm[u'orgo.arrowmodel'], null=False))
        ))
        db.create_unique(m2m_table_name, ['solutionmodel_id', 'arrowmodel_id'])

        # Adding model 'SynthesisProblemModel'
        db.create_table(u'orgo_synthesisproblemmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('solution', self.gf('django.db.models.fields.related.ForeignKey')(related_name='spSolution', null=True, to=orm['orgo.SolutionModel'])),
            ('target', self.gf('django.db.models.fields.related.ForeignKey')(related_name='spTarget', null=True, to=orm['orgo.MoleculeBoxModel'])),
            ('retain', self.gf('django.db.models.fields.BooleanField')()),
            ('solverCredited', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'orgo', ['SynthesisProblemModel'])

        # Adding M2M table for field molecules on 'SynthesisProblemModel'
        m2m_table_name = db.shorten_name(u'orgo_synthesisproblemmodel_molecules')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('synthesisproblemmodel', models.ForeignKey(orm[u'orgo.synthesisproblemmodel'], null=False)),
            ('moleculeboxmodel', models.ForeignKey(orm[u'orgo.moleculeboxmodel'], null=False))
        ))
        db.create_unique(m2m_table_name, ['synthesisproblemmodel_id', 'moleculeboxmodel_id'])

        # Adding M2M table for field arrows on 'SynthesisProblemModel'
        m2m_table_name = db.shorten_name(u'orgo_synthesisproblemmodel_arrows')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('synthesisproblemmodel', models.ForeignKey(orm[u'orgo.synthesisproblemmodel'], null=False)),
            ('arrowmodel', models.ForeignKey(orm[u'orgo.arrowmodel'], null=False))
        ))
        db.create_unique(m2m_table_name, ['synthesisproblemmodel_id', 'arrowmodel_id'])

        # Adding model 'ReactionStepModel'
        db.create_table(u'orgo_reactionstepmodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reactionStep', self.gf('orgo.models.PickledObjectField')()),
            ('reactantBox', self.gf('django.db.models.fields.related.ForeignKey')(related_name='reactant', null=True, on_delete=models.SET_NULL, to=orm['orgo.MoleculeBoxModel'])),
            ('productBox', self.gf('django.db.models.fields.related.ForeignKey')(related_name='product', null=True, on_delete=models.SET_NULL, to=orm['orgo.MoleculeBoxModel'])),
            ('html', self.gf('django.db.models.fields.TextField')(null=True)),
            ('done', self.gf('django.db.models.fields.BooleanField')()),
            ('catagory', self.gf('django.db.models.fields.CharField')(max_length=100, null=True)),
        ))
        db.send_create_signal(u'orgo', ['ReactionStepModel'])

        # Adding model 'ReagentType'
        db.create_table(u'orgo_reagenttype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
        ))
        db.send_create_signal(u'orgo', ['ReagentType'])

        # Adding model 'AccuracyModel'
        db.create_table(u'orgo_accuracymodel', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('catagory', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('correct', self.gf('django.db.models.fields.SmallIntegerField')()),
            ('total', self.gf('django.db.models.fields.SmallIntegerField')()),
        ))
        db.send_create_signal(u'orgo', ['AccuracyModel'])

        # Adding model 'UserProfile'
        db.create_table(u'orgo_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['auth.User'], unique=True)),
            ('currentNameReagentProblem', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['orgo.ReactionStepModel'], null=True, on_delete=models.SET_NULL)),
            ('currentSynthesisProblem', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['orgo.SynthesisProblemModel'], null=True, on_delete=models.SET_NULL)),
            ('autocompleteType', self.gf('django.db.models.fields.CharField')(max_length=20)),
            ('correctSynths', self.gf('django.db.models.fields.SmallIntegerField')(default=0)),
            ('assists', self.gf('django.db.models.fields.SmallIntegerField')(default=0)),
        ))
        db.send_create_signal(u'orgo', ['UserProfile'])

        # Adding M2M table for field savedReagentTypes on 'UserProfile'
        m2m_table_name = db.shorten_name(u'orgo_userprofile_savedReagentTypes')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('userprofile', models.ForeignKey(orm[u'orgo.userprofile'], null=False)),
            ('reagenttype', models.ForeignKey(orm[u'orgo.reagenttype'], null=False))
        ))
        db.create_unique(m2m_table_name, ['userprofile_id', 'reagenttype_id'])

        # Adding M2M table for field accuracies on 'UserProfile'
        m2m_table_name = db.shorten_name(u'orgo_userprofile_accuracies')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('userprofile', models.ForeignKey(orm[u'orgo.userprofile'], null=False)),
            ('accuracymodel', models.ForeignKey(orm[u'orgo.accuracymodel'], null=False))
        ))
        db.create_unique(m2m_table_name, ['userprofile_id', 'accuracymodel_id'])


    def backwards(self, orm):
        # Deleting model 'MoleculeBoxModel'
        db.delete_table(u'orgo_moleculeboxmodel')

        # Deleting model 'ArrowModel'
        db.delete_table(u'orgo_arrowmodel')

        # Deleting model 'SolutionModel'
        db.delete_table(u'orgo_solutionmodel')

        # Removing M2M table for field molecules on 'SolutionModel'
        db.delete_table(db.shorten_name(u'orgo_solutionmodel_molecules'))

        # Removing M2M table for field arrows on 'SolutionModel'
        db.delete_table(db.shorten_name(u'orgo_solutionmodel_arrows'))

        # Deleting model 'SynthesisProblemModel'
        db.delete_table(u'orgo_synthesisproblemmodel')

        # Removing M2M table for field molecules on 'SynthesisProblemModel'
        db.delete_table(db.shorten_name(u'orgo_synthesisproblemmodel_molecules'))

        # Removing M2M table for field arrows on 'SynthesisProblemModel'
        db.delete_table(db.shorten_name(u'orgo_synthesisproblemmodel_arrows'))

        # Deleting model 'ReactionStepModel'
        db.delete_table(u'orgo_reactionstepmodel')

        # Deleting model 'ReagentType'
        db.delete_table(u'orgo_reagenttype')

        # Deleting model 'AccuracyModel'
        db.delete_table(u'orgo_accuracymodel')

        # Deleting model 'UserProfile'
        db.delete_table(u'orgo_userprofile')

        # Removing M2M table for field savedReagentTypes on 'UserProfile'
        db.delete_table(db.shorten_name(u'orgo_userprofile_savedReagentTypes'))

        # Removing M2M table for field accuracies on 'UserProfile'
        db.delete_table(db.shorten_name(u'orgo_userprofile_accuracies'))


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'orgo.accuracymodel': {
            'Meta': {'object_name': 'AccuracyModel'},
            'catagory': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'correct': ('django.db.models.fields.SmallIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'total': ('django.db.models.fields.SmallIntegerField', [], {})
        },
        u'orgo.arrowmodel': {
            'Meta': {'object_name': 'ArrowModel'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'pointFrom': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'arrowPointsFrom'", 'null': 'True', 'to': u"orm['orgo.MoleculeBoxModel']"}),
            'pointTo': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'arrowPointsTo'", 'null': 'True', 'to': u"orm['orgo.MoleculeBoxModel']"}),
            'reagentsHtml': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        u'orgo.moleculeboxmodel': {
            'Meta': {'object_name': 'MoleculeBoxModel'},
            'equalsTarget': ('django.db.models.fields.BooleanField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'isStartingMaterial': ('django.db.models.fields.BooleanField', [], {}),
            'moleculeBox': ('orgo.models.PickledObjectField', [], {'null': 'True'}),
            'svg': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        u'orgo.reactionstepmodel': {
            'Meta': {'object_name': 'ReactionStepModel'},
            'catagory': ('django.db.models.fields.CharField', [], {'max_length': '100', 'null': 'True'}),
            'done': ('django.db.models.fields.BooleanField', [], {}),
            'html': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'productBox': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'product'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['orgo.MoleculeBoxModel']"}),
            'reactantBox': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'reactant'", 'null': 'True', 'on_delete': 'models.SET_NULL', 'to': u"orm['orgo.MoleculeBoxModel']"}),
            'reactionStep': ('orgo.models.PickledObjectField', [], {})
        },
        u'orgo.reagenttype': {
            'Meta': {'object_name': 'ReagentType'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'orgo.solutionmodel': {
            'Meta': {'object_name': 'SolutionModel'},
            'arrows': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'solutionArrows'", 'symmetrical': 'False', 'to': u"orm['orgo.ArrowModel']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecules': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'solutionMolecules'", 'symmetrical': 'False', 'to': u"orm['orgo.MoleculeBoxModel']"})
        },
        u'orgo.synthesisproblemmodel': {
            'Meta': {'object_name': 'SynthesisProblemModel'},
            'arrows': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'spArrows'", 'symmetrical': 'False', 'to': u"orm['orgo.ArrowModel']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molecules': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'spMolecules'", 'symmetrical': 'False', 'to': u"orm['orgo.MoleculeBoxModel']"}),
            'retain': ('django.db.models.fields.BooleanField', [], {}),
            'solution': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'spSolution'", 'null': 'True', 'to': u"orm['orgo.SolutionModel']"}),
            'solverCredited': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'target': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'spTarget'", 'null': 'True', 'to': u"orm['orgo.MoleculeBoxModel']"})
        },
        u'orgo.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            'accuracies': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['orgo.AccuracyModel']", 'symmetrical': 'False'}),
            'assists': ('django.db.models.fields.SmallIntegerField', [], {'default': '0'}),
            'autocompleteType': ('django.db.models.fields.CharField', [], {'max_length': '20'}),
            'correctSynths': ('django.db.models.fields.SmallIntegerField', [], {'default': '0'}),
            'currentNameReagentProblem': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['orgo.ReactionStepModel']", 'null': 'True', 'on_delete': 'models.SET_NULL'}),
            'currentSynthesisProblem': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['orgo.SynthesisProblemModel']", 'null': 'True', 'on_delete': 'models.SET_NULL'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'savedReagentTypes': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['orgo.ReagentType']", 'symmetrical': 'False'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        }
    }

    complete_apps = ['orgo']