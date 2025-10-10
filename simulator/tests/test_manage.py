from __future__ import annotations

from django.contrib.auth import get_user_model
from django.contrib.auth.models import Group
from django.test import TestCase
from django.urls import reverse

from clinic.models import Regimen
from simulator import models


class SimulatorManageViewsTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("viewer", password="pass123")
        self.editor = get_user_model().objects.create_user("editor", password="pass123")
        self.group = Group.objects.create(name="Simulator Editors")
        self.editor.groups.add(self.group)
        self.regimen_a = Regimen.objects.create(
            name="VRd",
            line="Frontline",
            components="Bortezomib, Lenalidomide, Dexamethasone",
            intent="curative",
            notes="Standard induction.",
        )
        self.regimen_b = Regimen.objects.create(
            name="Dara-Rd",
            line="Relapsed",
            components="Daratumumab, Lenalidomide, Dexamethasone",
            intent="salvage",
        )
        self.scenario = models.Scenario.objects.create(
            title="Newly diagnosed case",
            clinical_stage="newly_diagnosed",
            summary="Patient with symptomatic disease.",
            risk_stratification="Standard Risk",
            guideline_notes="Follow IMWG guidelines.",
            lab_snapshot={"IgG": "4500 mg/dL"},
        )
        self.scenario.recommended_regimens.add(self.regimen_a)

    def login_user(self, user) -> None:
        self.client.force_login(user)

    def test_regimen_create_requires_editor(self) -> None:
        self.login_user(self.user)
        url = reverse("simulator:regimen_create")
        response = self.client.post(
            url,
            {
                "name": "Ixazomib-Rd",
                "line": "Maintenance",
                "components": "Ixazomib, Lenalidomide, Dexamethasone",
                "intent": "maintenance",
                "notes": "Oral option.",
            },
        )
        self.assertEqual(response.status_code, 403)

    def test_editor_can_create_and_delete_regimen(self) -> None:
        self.login_user(self.editor)
        create_url = reverse("simulator:regimen_create")
        response = self.client.post(
            create_url,
            {
                "name": "KRd",
                "line": "Frontline",
                "components": "Carfilzomib, Lenalidomide, Dexamethasone",
                "intent": "curative",
                "notes": "",
            },
        )
        self.assertEqual(response.status_code, 302)
        regimen = Regimen.objects.get(name="KRd")

        delete_url = reverse("simulator:regimen_delete", args=[regimen.pk])
        response = self.client.get(
            delete_url,
            HTTP_HX_REQUEST="true",
        )
        self.assertContains(response, "Delete “KRd”", status_code=200)

        response = self.client.post(
            delete_url,
            {},
            HTTP_HX_REQUEST="true",
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn("HX-Redirect", response.headers)
        self.assertFalse(Regimen.objects.filter(name="KRd").exists())

    def test_editor_can_duplicate_scenario_with_regimens(self) -> None:
        self.scenario.recommended_regimens.add(self.regimen_b)
        self.login_user(self.editor)
        duplicate_url = reverse("simulator:scenario_duplicate", args=[self.scenario.pk])
        response = self.client.post(duplicate_url)
        self.assertEqual(response.status_code, 302)
        duplicate = models.Scenario.objects.get(title=f"{self.scenario.title} (copy)")
        self.assertCountEqual(
            duplicate.recommended_regimens.values_list("pk", flat=True),
            self.scenario.recommended_regimens.values_list("pk", flat=True),
        )

    def test_htmx_add_and_remove_regimen(self) -> None:
        self.login_user(self.editor)
        add_url = reverse("simulator:scenario_regimen_add", args=[self.scenario.pk])
        response = self.client.post(
            add_url,
            {"regimen_id": [self.regimen_b.pk]},
            HTTP_HX_REQUEST="true",
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn("hx-swap-oob", response.content.decode())
        self.assertTrue(self.scenario.recommended_regimens.filter(pk=self.regimen_b.pk).exists())

        remove_url = reverse("simulator:scenario_regimen_remove", args=[self.scenario.pk])
        response = self.client.post(
            remove_url,
            {"regimen_id": self.regimen_a.pk},
            HTTP_HX_REQUEST="true",
        )
        self.assertEqual(response.status_code, 200)
        self.assertFalse(self.scenario.recommended_regimens.filter(pk=self.regimen_a.pk).exists())

    def test_non_editor_htmx_add_denied(self) -> None:
        self.login_user(self.user)
        add_url = reverse("simulator:scenario_regimen_add", args=[self.scenario.pk])
        response = self.client.post(
            add_url,
            {"regimen_id": [self.regimen_b.pk]},
            HTTP_HX_REQUEST="true",
        )
        self.assertEqual(response.status_code, 403)
